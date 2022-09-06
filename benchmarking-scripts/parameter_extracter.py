import sys
import subprocess
import argparse
import statistics
import os
import glob
import shutil
import math
# from scipy.stats import shapiro
# import matplotlib.pyplot as plt

def process_params(clamr_output, proc_num, exe):
    """
    clamr_output: output from clamr executable which is being processed
    proc_num:     number of processes used for the mpiexec launch
    exe:          name of the clamr executable being run
    """

    nowned = []
    nremote = []
    comm_partners = []
    indices_needed = []
    nsizes = []
    block_sizes = []

    for line in output.splitlines():
        if "nowned" in line:
            line = line.split("-")
            if line[1].strip().isdigit():
                nowned.append(int(line[1]))
        elif "nremote" in line:
            line = line.split("-")
            if line[1].strip().isdigit():
                nremote.append(int(line[1]))
        elif "num_comm_partners" in line:
            line = line.split("-")
            if line[1].strip().isdigit():
                comm_partners.append(int(line[1]))
        elif "indices_needed" in line:
            line = line.split("-")[1].strip().split(" ")

            # initialize nested lists of indices_needed
            # if indices_needed is empty
            if len(indices_needed) == 0:
                indices_needed = [[]] * len(line)

            # append value to nested list
            for i, val in enumerate(line):
                indices_needed[i].append(int(line[i]))
        elif "nsizes" in line:
            line = line.split("-")
            line = line[1].strip().split(" ")

            # initialize nested lists of indices_needed
            # if indices_needed is empty
            if len(nsizes) == 0:
                nsizes = [[]] * len(line)

            # append value to nested list
            for i, val in enumerate(line):
                nsizes[i].append(int(line[i]))
        elif "blocksize" in line:
            line = line.split("-")
            if line[1].strip().isdigit():
                block_sizes.append(int(line[1]))
        elif "PARAM:" in line:
            print(line)


    # create directory for each process count graphs
    dir_path = "./param_info_" + str(proc_num) + "/"
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    print("INFO: for the following parameters, a p-value of > 0.05 implies a gaussian distribution\n")

    print("nowned: " + str(round(statistics.mean(nowned))))
    print("nowned stdev: " + str(round(statistics.stdev(nowned))))
    #print("nowned p-value: " + str(shapiro(nowned).pvalue) + "\n")
    # plt.hist(nowned, bins=20)
    # plt.title("nowned size (" + str(proc_num) + " processes)")
    # plt.xlabel("size (bytes)")
    # plt.ylabel("frequency")
    # plt.savefig(dir_path + "/" + exe + "_nowned.png")
    # plt.clf()

    print("nremote: " + str(round(statistics.mean(nremote))))
    print("nremote stdev: " + str(round(statistics.stdev(nremote))))
    #print("nremote p-value: " + str(shapiro(nremote).pvalue) + "\n")
    # plt.hist(nremote, bins=20)
    # plt.title("nremote size (" + str(proc_num) + " processes)")
    # plt.xlabel("size (bytes)")
    # plt.ylabel("frequency")
    # plt.savefig(dir_path + "/" + exe + "_nremote.png")
    # plt.clf()

    print("num_comm_partners: " + str(round(statistics.mean(comm_partners))) + "\n")
    # plt.hist(comm_partners, bins=20)
    # plt.title("comm_partners size (" + str(proc_num) + " processes)")
    # plt.xlabel("number of partners")
    # plt.ylabel("frequency")
    # plt.savefig(dir_path + "/" + exe + "_comm_partners.png")
    # plt.clf()

    if len(block_sizes) == 0:
        print("Error: Blocksize was not printed for proccess count " + str(proc_num) + "\n")
    else:
        print("block_sizes: " + str(round(statistics.mean(block_sizes))))
        print("block_sizes stdev: " + str(round(statistics.stdev(block_sizes))))
        #print("block_sizes p-value: " + str(shapiro(block_sizes).pvalue))
        # plt.hist(block_sizes, bins=20)
        # plt.title("block_sizes size (" + str(proc_num) + " processes)")
        # plt.xlabel("size (bytes)")
        # plt.ylabel("frequency")
        # plt.savefig(dir_path + "/" + exe + "_block_sizes.png")
        # plt.clf()

    return [
        round(statistics.mean(nowned)),
        round(statistics.stdev(nowned)),
        round(statistics.mean(nremote)),
        round(statistics.stdev(nremote)),
        round(statistics.mean(comm_partners)),
        round(statistics.mean(block_sizes)),
        round(statistics.stdev(block_sizes))
    ]

    # pvalue_list = []
    # for index, val in enumerate(indices_needed):
    #     #pvalue_list.append(shapiro(indices_needed[index]).pvalue)
    #     indices_needed[index] = round(statistics.mean(indices_needed[index]))
    # print("indices_needed: " + str(indices_needed))
    # #print("indices_needed p-value: " + str(pvalue_list) + "\n")

    # pvalue_list = []
    # for index, val in enumerate(nsizes):
    #     #pvalue_list.append(shapiro(nsizes[index]).pvalue)
    #     nsizes[index] = round(statistics.mean(nsizes[index]))
    # print("nsizes: " + str(nsizes) + "\n\n")
    # #print("nsizes p-value: " + str(pvalue_list) + "\n\n")

def run_clamr_with_params(proc_num, nowned_avg, nowned_stdv,
                          nremote_avg, nremote_stdv, comm_partners,
                          block_size_avg=3, block_size_stdv=3):
    # If process count is greater than logical processes
    if nproc < proc_num:
        cmd = ['mpirun',
         '--oversubscribe',
         '--use-hwthread-cpus',
         '-np',
         str(proc_num),
         './clamr/build/l7_update_perf',
         '-o', str(nowned_avg),
         '-O', str(nowned_stdv),
         '-r', str(nremote_avg),
         '-R', str(nremote_stdv),
         '-b', str(block_size_avg),
         '-B', str(block_size_stdv),
         '-n', str(comm_partners),
         '-I', str(10),
     ]
    else:
        cmd = ['mpirun', '-np', str(proc_num), './clamr/build/l7_update_perf',
               '-o', str(nowned_avg), '-r', str(nremote_avg),
               '-n', str(comm_partners), '-b', str(block_size_avg), '-I', str(10)]

    output = subprocess.run(cmd,
                            capture_output=False, encoding='UTF-8').stdout

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_process',
                        dest='proc_num',
                        default=-1,
                        action='store',
                        nargs='?',
                        type=int,
                        help='Specify the maximum number of processes used')

    parser.add_argument('-c', '--clean',
                        action='store_true',
                        help='Removes previously generated files')

    parser.add_argument('-r', '--run',
                        action='store_true',
                        help='Run even when clean is active')

    args = parser.parse_args()
    if sys.platform == 'darwin':
        nproc = int(subprocess.run(['sysctl', '-n', 'hw.ncpu'], capture_output=True, encoding='UTF-8').stdout.strip())
    elif sys.platform == 'linux':
        nproc = int(subprocess.run(['nproc'], capture_output=True, encoding='UTF-8').stdout.strip())
    else:
        nproc = 4

    if args.clean:
        file_list = glob.glob('./out*', recursive=True)
        file_list.extend(glob.glob('./param_info_*/', recursive=True))
        for file_path in file_list:
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
                else:
                    print("Error unknown file type at : ", file_path)
            except:
                print("Error while deleting file : ", file_path)
        if os.path.isdir('./clamr'):
            shutil.rmtree('./clamr')
        if not args.run:
            exit()

    # if process count is unset, use max process count
    # otherwise, set from argparse
    if args.proc_num == -1:
        proc_num = nproc
    else:
        proc_num = args.proc_num

    # bootstrap clamr (with code to extract parameters)
    clamr_path = "./clamr/"
    if not os.path.exists(clamr_path):
        print("Boostrapping CLAMR")
        subprocess.run(['git', 'clone', '-b', 'parameters',
                        'https://github.com/carsonwoods/clamr'],
                       capture_output=True, encoding='UTF-8')

        cwd = os.getcwd() # remember our original working directory
        os.mkdir('./clamr/build')
        os.chdir('./clamr/build/')
        subprocess.run(['cmake', '..'],
           capture_output=True, encoding='UTF-8')
        subprocess.run(['make', '-j', str(nproc)],
           capture_output=True, encoding='UTF-8')
        os.chdir(cwd) # get back to our original working directory
        print("Done")


    #clamr_name = ['clamr_mpionly', 'clamr_mpicheck', 'clamr_mpiopenmponly']
    clamr_name = ['clamr_mpionly', 'clamr_mpicheck']

    for exe in clamr_name:
        print_str = " Running " + exe + " on " + str(proc_num) + " processes "
        target_size = os.get_terminal_size().columns
        diff_size = target_size - len(print_str)
        if diff_size <= 0:
            print(print_str)
        else:
            for x in range(0,round(diff_size/2)):
                if target_size-len(print_str) == 1:
                    print_str += "#"
                    break
                elif target_size-len(print_str) == -1:
                    print_str = print_str[:len(print_str)-1]
                    break
                elif target_size == len(print_str):
                    break
                else:
                    print_str = "#" + print_str + "#"

            print(print_str)

        # If process count is greater than logical processes
        if nproc < proc_num:
            clamr_cmd = ['mpirun',
             '--oversubscribe',
             '--use-hwthread-cpus',
             '-np',
             str(proc_num),
             "./clamr/build/" + exe]
        else:
            clamr_cmd = ['mpirun', '-np', str(proc_num), "./clamr/build/"+ exe]

        output = subprocess.run(clamr_cmd,
                                capture_output=True, encoding='UTF-8').stdout
        params = process_params(output, proc_num, exe)

        # blocksize is not always reported
        # this ensures that, if not reported, parameters are not passed
        if len(params) > 4:
            run_clamr_with_params(proc_num=proc_num, nowned_avg=params[0],
                                  nowned_stdv=params[1], nremote_avg=params[2],
                                  nremote_stdv=params[3], comm_partners=params[4],
                                  block_size_avg=params[5], block_size_stdv=params[5])
        else:
            run_clamr_with_params(proc_num=proc_num, nowned_avg=params[0],
                                  nowned_stdv=params[1], nremote_avg=params[2],
                                  nremote_stdv=params[3], comm_partners=params[4])


    file_list = glob.glob('./out*', recursive=True)
    for file_path in file_list:
        try:
            os.remove(file_path)
        except:
            print("Error while deleting file : ", file_path)
