#!/usr/bin/python3.6
import sys
import os

cwd = os.getcwd()

input_file = open(sys.argv[1], 'r')
Lines = input_file.readlines()

count = 0
Dict = {}

Dict["Isend_max"] = [0.0,0] 
Dict["Isend_min"] = [0.0,0]  
Dict["Isend_avg"] = [0.0,0] 

Dict["Irecv_max"] = [0.0,0] 
Dict["Irecv_min"] = [0.0,0] 
Dict["Irecv_avg"] = [0.0,0] 

Dict["Num_partners_max"] = [0.0,0] 
Dict["Num_partners_min"] = [0.0,0] 
Dict["Num_partners_avg"] = [0.0,0] 


for line in Lines:
  count += 1

  if "Isend max:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Isend_max"][0]+=float(value)
    Dict["Isend_max"][1]+=1


  elif "Isend min:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Isend_min"][0]+=float(value)
    Dict["Isend_min"][1]+=1


  elif "Isend avg:" in line :
    start_index = line.find(":")+1
    value = line[start_index : ].strip() 
    Dict["Isend_avg"][0]+=float(value)
    Dict["Isend_avg"][1]+=1


  elif "Irecv max:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Irecv_max"][0]+=float(value)
    Dict["Irecv_max"][1]+=1


  elif "Irecv min:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Irecv_min"][0]+=float(value)
    Dict["Irecv_min"][1]+=1


  elif "Irecv avg:" in line :
    start_index = line.find(":")+1
    value = line[start_index : ].strip()
    Dict["Irecv_avg"][0]+=float(value)
    Dict["Irecv_avg"][1]+=1

  elif "Num partners max:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Num_partners_max"][0]+=int(value)
    Dict["Num_partners_max"][1]+=1

  elif "Num partners min:" in line :
    start_index = line.find(":")+1
    end_index = line.find("(")-1
    value = line[start_index : end_index].strip()
    Dict["Num_partners_min"][0]+=int(value)
    Dict["Num_partners_min"][1]+=1

  elif "Num partners avg:" in line :
    start_index = line.find(":")+1
    value = line[start_index : ].strip()
    Dict["Num_partners_avg"][0]+=float(value)
    Dict["Num_partners_avg"][1]+=1




Isend_max_avg = float(Dict["Isend_max"][0] / Dict["Isend_max"][1])
Isend_min_avg = float(Dict["Isend_min"][0] / Dict["Isend_min"][1])
Isend_avg = float(Dict["Isend_avg"][0] / Dict["Isend_avg"][1])
Irecv_max_avg = float(Dict["Irecv_max"][0] / Dict["Irecv_max"][1])
Irecv_min_avg = float(Dict["Irecv_min"][0] / Dict["Irecv_min"][1])
Irecv_avg = float(Dict["Irecv_avg"][0] / Dict["Irecv_avg"][1])
Num_partners_max_avg = float(Dict["Num_partners_max"][0] / float(Dict["Num_partners_max"][1]))
Num_partners_min_avg = float(Dict["Num_partners_min"][0] / float(Dict["Num_partners_min"][1]))
Num_partners_avg = float(Dict["Num_partners_avg"][0] / float(Dict["Num_partners_avg"][1]))



print("\n",cwd);
print("\nIsend max avg: {:.10f}\n".format(Isend_max_avg))
print("\nIsend min avg: {:.10f}\n".format(Isend_min_avg))
print("\nIsend avg: {:.10f}\n".format(Isend_avg))
print("\nIrecv max avg: {:.10f}\n".format(Irecv_max_avg))
print("\nIrecv min avg: {:.10f}\n".format(Irecv_min_avg))
print("\nIrecv avg: {:.10f}\n".format(Irecv_avg))
print("\nNum partners max avg: {}\n".format(Num_partners_max_avg))
print("\nNum partners min avg: {}\n".format(Num_partners_min_avg))
print("\nNum partners avg: {}\n".format(Num_partners_avg))
