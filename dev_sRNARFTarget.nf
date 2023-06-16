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

  chunksize = 500
  i=0
  dictionary_list=[]
  for chunk in pd.read_csv("$srna", engine='c', chunksize=chunksize, skiprows=lambda x: (x != 0) and not x % 2):
    print("Completed", chunksize * i, flush=True)
    for index, row in chunk.iterrows():
        s = Sequence(row[0])
        freqs = s.kmer_frequencies(3, relative=True, overlap=True)
        dictionary_list.append(freqs)
    i=i+1
  
  # Deprecated
  #df = pd.DataFrame(columns=kmer_combinations)
  #df = df.append(dictionary_list,ignore_index=True).fillna(0)
  
  df = pd.concat([pd.DataFrame(columns=kmer_combinations), pd.DataFrame.from_dict(dictionary_list)]).fillna(0)
  df = pd.DataFrame(np.repeat(df.values, count, axis=0), columns=kmer_combinations)
  df = df.round(9) # We round to 5 in last process. Rounding now reduces size and increase performance
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
  df = pd.DataFrame(columns=kmer_combinations)
  chunksize = 100
  seqarr = []
  i=0
  for chunk in pd.read_csv("$mrna", chunksize=chunksize, header=0, skiprows=lambda x: (x != 0) and not x % 2):
      print("Completed", chunksize * i, flush=True)
      for index, row in chunk.iterrows():
          s = Sequence(row[0])
          freqs = s.kmer_frequencies(3, relative=True, overlap=True)
          seqarr.append(freqs)
          df = pd.concat([df, pd.DataFrame(freqs, index=[0])], ignore_index=True).fillna(0)
          df = df.round(9)
      i=i+1
  
  header = True
  for x in range(scount):
      if x % 100 == 0:
          print(x , "of" , scount, flush=True)
      df.to_csv(fileout, header=header, mode='a', index=False)
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
  file diffpkl from setResult3B
  file lrf from ldedmodel

  output:
  file 'Results_pred_probs.npy' into process5result

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import numpy  as np
  import pickle
  import os
  import gc

  fileout="Results_pred_probs"
  if os.path.exists(fileout + ".npy"):
    os.remove(fileout + ".npy")

  diff = pd.read_pickle(diffpkl)

  # load the saved random forest model from disk
  loaded_RFmodel = pickle.load(open(lrf, 'rb'))

  #predict probabilities for class 1
  predict_proba = loaded_RFmodel.predict_proba(diff.fillna(0))

  del loaded_RFmodel
  del testdf 
  gc.collect()

  #write probabilities to file
  np.save(fileout, predict_proba[:, 1], allow_pickle=True, fix_imports=True)
  """
}
process5result.set{setResult5}

//-------------------------Process_5---------------------------//

process generateSortedResultFile{

  input:
  file mlfile from setResult5
  file ns3file from setResult111
  file difffile from setResult33B

  output:
  file 'Prediction_probabilities.csv' into process6result1
  file 'FeatureFile.csv' into process6result2

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import pickle
  import gc

  #Generate sorted prediction result file
  df1 = pd.read_pickle("$ns3file").iloc[:, 0:2] 
  df2 = pd.DataFrame(np.load(new, mmap_mode=None, allow_pickle=True)).round(5)
  
  #TODO: As we round by 5, can we change dtype from default float64 to float32
  #df2 = pd.DataFrame(np.load(new, mmap_mode=None, allow_pickle=True)).round(5).astype('float32')

  df3 = pd.DataFrame(data=df1.values,columns=['sRNA_ID', 'mRNA_ID']).assign(Prediction_Probability=df2)

  del df2

  df4 = df3.sort_values('Prediction_Probability',ascending=False)
  del df3

  df4.to_csv('Prediction_probabilities.csv', sep='\t', index=False)
  del df4
  gc.collect()

  #Generate feature file with pair ids to be used for interpretability later
  dfp63 = pd.concat([df1, pd.read_pickle("$difffile")], axis = 1)
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
