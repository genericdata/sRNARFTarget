#!/usr/bin/env nextflow

//-------------------------channels----------------------------//
params.s = "srna.fasta"
params.m = "mrna.fasta"

srna = file("${params.s}")
mrna = file("${params.m}")

ldedmodel = Channel.fromPath('./PickledModelData/RFModel/sRNARFTargetModel.pickle')
if(!srna.exists()) {
exit 1, "The specified input sRNA fasta file does not exist: ${params.s}"
}

if(!mrna.exists()) {
exit 1, "The specified input mRNA fasta file does not exist: ${params.m}"
}

//------------------------Process_1---------------------------//

process createAllPossiblePairs{

  output:
  file 'pairs_names_seqs.pkl' into process1result mode flatten

  script:
  """
  #!/usr/bin/env python3
  import pickle
  import pandas as pd
  from Bio import SeqIO
  import os
  import pyarrow

  fileout="pairs_names_seqs.pkl"
  if os.path.exists(fileout):
    os.remove(fileout)

  records_s = list(SeqIO.parse("$srna", "fasta"))
  records_m = list(SeqIO.parse("$mrna", "fasta"))

  df = pd.DataFrame(columns=["sRNA_ID","mRNA_ID","sRNA_Sequence","mRNA_Sequence"])
  dictionary_list = []

  for i in range(len(records_s)):
      for j in range(len(records_m)):
          dictionary_data = {"sRNA_ID": records_s[i].id, "mRNA_ID": records_m[j].id, "sRNA_Sequence": str(records_s[i].seq), "mRNA_Sequence": str(records_m[j].seq)}
          dictionary_list.append(dictionary_data)

  df = pd.DataFrame.from_dict(dictionary_list)
  df.to_pickle(fileout)

  """
}
process1result.into{setResult1; setResult11; setResult111}

//-------------------------Process_2---------------------------//

process getsRNATrinucleotidesFrequncies{

  output:
  file 'sRNA_3mer.pkl' into process2result

  script:
  """
  #!/usr/bin/env python3
  from skbio import Sequence
  from itertools import product
  import pandas as pd
  import numpy as np
  import pickle
  import os

  def versiontuple(v):
      return tuple(map(int, (v.split("."))))

  def all_kmer_subsets(ss=["A", "T", "G", "C"]):
      return [''.join(p) for p in product(ss, repeat=3)]

  kmer_combinations = all_kmer_subsets()

  with open("$mrna", 'r') as fp:
      for count, line in enumerate(fp):
          pass
  count=int((count+1)/2)

  fileout="sRNA_3mer.pkl"
  if os.path.exists(fileout):
    os.remove(fileout)

  df = pd.DataFrame(columns=kmer_combinations)

  if versiontuple(pd.__version__) >= versiontuple("1.4.0"):
     engine="pyarrow"
  else:
     engine="c"
     
  chunksize = 500
  i=0
  dictionary_list=[]
  for chunk in pd.read_csv("$srna", engine=engine, chunksize=chunksize, skiprows=lambda x: (x != 0) and not x % 2):
    print("Completed", chunksize * i, flush=True)
    for index, row in chunk.iterrows():
        s = Sequence(row[0])
        freqs = s.kmer_frequencies(3, relative=True, overlap=True)
        dictionary_list.append(freqs)
    i=i+1

  df = df.append(dictionary_list,ignore_index=True).fillna(0)
  df = pd.DataFrame(np.repeat(df.values, count, axis=0), columns=kmer_combinations)
  df = df.round(8) # We round to 5 in last process. Rounding now reduces size and increase performance
  df.to_pickle(fileout)

  """
}

//Collect file
process2result.set{setResult2}

//-------------------------Process_3---------------------------//
process getmRNATrinucleotidesFrequncies{

  output:
  file 'mRNA_3mer.txt' into process3result

  script:
  """
  #!/usr/bin/env python3
  from skbio import Sequence
  from itertools import product
  import pandas as pd
  #import numpy as np

  def all_kmer_subsets(ss=["A", "T", "G", "C"]):
      return [''.join(p) for p in product(ss, repeat=3)]

  kmer_combinations = all_kmer_subsets()

  with open("$srna", 'r') as fp:
      for scount, line in enumerate(fp):
          pass
  scount=int((scount+1)/2)

  with open("$mrna", 'r') as fp:
      for mcount, line in enumerate(fp):
          pass
  mcount=(mcount+1)/2

  fileout="mRNA_3mer.txt"
  # Create a 0x64 dataframe
  dfn = pd.DataFrame(columns=kmer_combinations)
  header = True
  chunksize = 100
  seqarr = []
  i=0
  for chunk in pd.read_csv("$mrna", chunksize=chunksize, header=0, skiprows=lambda x: (x != 0) and not x % 2):
      print("Completed", chunksize * i, flush=True)
      for index, row in chunk.iterrows():
          s = Sequence(row[0])
          freqs = s.kmer_frequencies(3, relative=True, overlap=True)
          seqarr.append(freqs)
          dfn = pd.concat([dfn, pd.DataFrame(freqs, index=[0])], ignore_index=True).fillna(0)
      i=i+1

  print("dataframe =",dfn.shape, flush=True)
  print("Writing dataframe to", fileout, scount, "times", flush=True)

  for x in range(scount):
      if x % 100 == 0:
          print(x , "of" , scount, flush=True)
      dfn.to_csv(fileout, header=header, mode='a', index=False)
      header = False
"""
}
//Collect file
process3result.into{setResult3; setResult33}

//-------------------------Process_3B---------------------------//
process getDifference{

  input:
  file sRNas from setResult2
  file mRNas from setResult3

  output:
  file '3merdifference.pkl' into process3Bresult

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import pyarrow
  import gc
  import os
  import pickle

  fileout="3merdifference.pkl"
  if os.path.exists(fileout):
    os.remove(fileout)
    
  def versiontuple(v):
    return tuple(map(int, (v.split("."))))

  # Read in mRNAdf
  if versiontuple(pd.__version__) >= versiontuple("1.4.0"):
    engine="pyarrow"
  else:
    engine="c"
  
  mRNAdf = pd.read_csv(mRNas, engine=engine)
  
  sRNAdf = pd.read_pickle(sRNas)
  diff = mRNAdf.subtract(sRNAdf)

  # Clean up
  del mRNAdf
  del sRNAdf
  gc.collect()

  # Write out to pickle
  diff.to_pickle(fileout)
  #output8.to_csv(fileout, header=True, index=False, sep='\t', mode='a')

"""

}
//Collect file
process3Bresult.into{setResult3B; setResult33B}

//-------------------------Process_4---------------------------//

process runRandomForestModel{

input:
file newdata from setResult3
file lrf from ldedmodel

output:
file 'Results_pred_probs.txt' into process5result

script:
"""
#!/usr/bin/env python3
import pandas as pd
import numpy  as np
import pickle

def pred_prob(testdf):
    # load the saved random forest model from disk
    loaded_RFmodel = pickle.load(open('$lrf', 'rb'))
    #predict probabilities for class 1
    predict_proba = loaded_RFmodel.predict_proba(testdf)
    print(predict_proba)
    #write probabilities to file
    for i in predict_proba:
        with open('Results_pred_probs.txt','a') as fd:
            fd.write(str(i[1])+"\\n")
    print(predict_proba[:, 1])
    return predict_proba[:, 1]

testdata = pd.read_csv('$newdata', sep='\t', header=0)
testdf = pd.DataFrame(data = testdata)
testdf = testdf.fillna(0)

pred_prob(testdf)
"""
}
process5result.set{setResult5}

//-------------------------Process_5---------------------------//

process generateSortedResultFile{

input:
file mlfile from setResult5
file ns3file from setResult111
file difffile from setResult33

output:
file 'Prediction_probabilities.csv' into process6result1
file 'FeatureFile.csv' into process6result2

script:
"""
#!/usr/bin/env python3
import pandas as pd

#Generate sorted prediction result file
df1 = pd.read_csv('$ns3file', sep='\t', header=0)
df2 = pd.read_csv('$mlfile', sep='\t', header=None)
df3 = pd.DataFrame(data=df1.iloc[:, 0:2].values,columns=['sRNA_ID', 'mRNA_ID']).assign(Prediction_Probability=df2.round(5))

df4 = df3.sort_values('Prediction_Probability',ascending=False)
df4.to_csv('Prediction_probabilities.csv', sep='\t', index=False)

#Generate feature file with pair ids to be used for interpretability later
dfp61 = pd.read_csv('$difffile', sep='\t',header = 0) 
dfp63 = pd.concat([df1.iloc[:, 0:2], dfp61], axis = 1)
dfp63.to_csv('FeatureFile.csv', header = True, sep='\t', index=False)

"""
}

process6result1.collectFile(name: 'Prediction_probabilities.csv', storeDir:'sRNARFTargetResult')
process6result2.collectFile(name: 'FeatureFile.csv', storeDir:'sRNARFTargetResult')


//-------------------------summary---------------------------//

workflow.onComplete {
println(
"""
Pipeline execution summary
---------------------------
Run as : ${workflow.commandLine}
Completed at: ${workflow.complete}
Duration : ${workflow.duration}
Success : ${workflow.success}
workDir : ${workflow.workDir}
exit status : ${workflow.exitStatus}
""")
}
