import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(400,500,502,504),
    session=None):
    session=session or requests.Session()
    retry=Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,

    )
    adapter=HTTPAdapter(max_retries=retry)
    session.mount('http://',adapter)
    session.mount('https://',adapter)
    return session

def get_genesequence():
    import time
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/"+ ensembl_chr +":"+ ensembl_start + ".." + ensembl_end + "?mask=soft"  #1:55505221..55507220?mask=soft"
    t0=time.time()
    try:
        sequence = requests_retry_session().get(server+ext, headers={"Content-Type" : "text/x-fasta"})
    except Exception as x:
        print('It failed :(', x._class_._name_)
    else:
        print('It eventually worked', sequence.status_code)
    finally:
        t1=time.time()
        print('Took',t1-t0,'seconds')

    fasta=sequence.text.upper()
    with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/"+ fasta_file,"w") as fastafile:
        fastafile.write(fasta)




def selection():
    import subprocess
    import os

    commands = "python dsNickFury3.3.py -m selection -s 20_NGG -g hg38 --cluster --matchSiteCutoff 5000 --targetFasta " + fasta_file + " --outputToFile "+ "/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/"+selection_output
    proc=subprocess.call(commands,shell=True)
    print(proc)




# write header
def result_file():
    import csv
    with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/"+selection_result , "w") as result:
        header = ['Guide', 'Mismatch Risk', 'Azimuth Score', 'Chr', 'Start', 'End', 'Gene', 'Strand', 'Mismatches: 0',
                  'Mismatches: 1', 'Mismatches: 2', 'Mismatches: 3', "Off target"]
        writer = csv.writer(result, delimiter="\t")
        writer.writerow(header)





# divide the result per guide
def result_chunk():
    import re

    with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/"+ selection_output, "r") as output_file:
        line_output = [line.strip().split('\n') for line in output_file]
        info_list = []
        decisioner = 1
        info_str = []
        line_output.append(['GATTAATATCCTTGTCATAA_AGG\tMismatch Risk: 294'])
        for line in line_output:
            guide_regexmatch = re.search('(\S{24})\tMismatch', line[0])
            if (guide_regexmatch):
                decisioner = decisioner * (-1)
            if (decisioner == -1):
                info_str.append(line[0])
            if (decisioner == 1):
                info_list.append(info_str)
                info_str = [line[0]]
                decisioner = -1
    info_list.pop(0)
    info_list.pop(0)
    return (info_list)


def selection__result():
    import re
    import csv


    info_matrix = result_chunk()
    result_dictionary = {'Guide': '', 'Mismatch Risk': '', 'Azimuth Score': '', 'Chr': '', 'Start': '', 'End': '', 'Gene':'', 'Strand': '','Mismatches: 0': '', 'Mismatches: 1': '',
                         'Mismatches: 2': '', 'Mismatches: 3': '', 'Off target': ''}

    for i in range(0, len(info_matrix)):  # len(info_matrix)
        count_mismatch0 = 0
        count_mismatch1 = 0
        count_mismatch2 = 0
        count_mismatch3 = 0
        decisioner0 = -1
        decisioner1 = -1
        decisioner2 = -1
        decisioner3 = -1
        info_list = info_matrix[i]
        for j in range(0, len(info_list)):

            breakpoint_regexmatch = re.search("(\S{24})\tMismatch Risk: 9999999", info_list[j])
            if breakpoint_regexmatch:
                break
            # print(info_list[j])

            # match for guide
            guide_regexmatch = re.search('(\S{24})\tMismatch', info_list[j])
            if (guide_regexmatch):
                result_dictionary["Guide"] = guide_regexmatch.group(1)

            # match for risk
            risk_regexmatch = re.search('Risk:\s(.+)', info_list[j])
            if (risk_regexmatch):
                result_dictionary["Mismatch Risk"] = risk_regexmatch.group(1)
            elif (result_dictionary["Mismatch Risk"] == ''):
                result_dictionary["Mismatch Risk"] = 'None'

            # match for azimuth score
            azimuth_regexmatch = re.search("Score:\s(\d.\d*)", info_list[j])
            if (azimuth_regexmatch):
                result_dictionary["Azimuth Score"] = azimuth_regexmatch.group(1)


            elif (result_dictionary["Azimuth Score"] == ''):
                result_dictionary["Azimuth Score"] = 'None'

            # count the number of "Mismatches :1"
            if info_list[j] == "Mismatches: 1":
                decisioner1 = decisioner1 * -1
            if info_list[j] == "Mismatches: 2":
                decisioner1 = decisioner1 * -1
            if decisioner1 == 1:
                count_mismatch1 = count_mismatch1 + 1
            result_dictionary["Mismatches: 1"] = count_mismatch1 - 1

            # count the number of "Mismatches :2"
            if info_list[j] == "Mismatches: 2":
                decisioner2 = decisioner2 * -1
            if info_list[j] == "Mismatches: 3":
                decisioner2 = decisioner2 * -1
            if decisioner2 == 1:
                count_mismatch2 = count_mismatch2 + 1
            result_dictionary["Mismatches: 2"] = count_mismatch2 - 1

            # count the number of "Mismatches: 3"
            if info_list[j] == "Mismatches: 3":
                decisioner3 = decisioner3 * -1
            if decisioner3 == 1:
                count_mismatch3 = count_mismatch3 + 1
            result_dictionary["Mismatches: 3"] = count_mismatch3 - 1

            # count the number of "Mismatches: 0"
            if info_list[j] == "Mismatches: 0":
                decisioner0 = decisioner0 * -1
            if info_list[j] == "Mismatches: 1":
                decisioner0 = decisioner0 * -1
            if decisioner0 == 1:
                count_mismatch0 = count_mismatch0 + 1


                #Check if the guide have perfect off target
                chr_regexmatch = re.search('(.+?)\t.+\t.+\t.+\t.+\t[+|-]', info_list[j])
                start_regexmatch = re.search(".+\t(.+?)\t.+\t.+\t.+\t[+|-]", info_list[j])
                end_regexmatch = re.search(".+\t.+\t(.+?)\t.+\t.+\t[+|-]", info_list[j])
                guide_check_regexmatch = re.search(".+\t.+\t.+\t.+(\S{24})\t.+\t[+|-]", info_list[j]) #To check which guide is the guide that in the title of the chunk, added because found out that when there are perfect off targets and the off target are in the same gene, the parser would record the last start and end position previously. So add this line to check the guide
                gene_regexmatch = False
                strand_regexmatch = False
                if chr_regexmatch and start_regexmatch and end_regexmatch and guide_check_regexmatch:
                    chr_result = chr_regexmatch.group(1)
                    start_result = start_regexmatch.group(1)
                    end_result = end_regexmatch.group(1)
                    guide_check = guide_check_regexmatch.group(1)
                    if chr_result == ensembl_chr and int(start_result) >= int(ensembl_start) and int(end_result)<= int(ensembl_end) and guide_check == guide_regexmatch.group(1):
                        result_dictionary["Chr"] = chr_result
                        result_dictionary["Start"] = start_result
                        result_dictionary["End"] = int(end_result)-2
                        result_dictionary["Mismatches: 0"] = count_mismatch0 - 2
                        gene_regexmatch = re.search(".+\t.+\t.+\t(.+)/", info_list[j])
                        strand_regexmatch = re.search(".+\t.+\t.+\t.+\t.+\t(.)\t\t\t", info_list[j])
                        if count_mismatch0-2 == 0:
                            result_dictionary["Off target"] = "No"
                        else:
                            result_dictionary["Off target"] = "Yes"

                    else:
                        gene_regexmatch=False
                        strand_regexmatch=False
                        result_dictionary["Mismatches: 0"] = count_mismatch0 - 1
                        if count_mismatch0-1 == 0:
                            result_dictionary["Off target"] = "No"
                        else:
                            result_dictionary["Off target"] = "Yes"




                if gene_regexmatch:
                    result_dictionary["Gene"] = gene_regexmatch.group(1)



                if strand_regexmatch:
                    result_dictionary["Strand"]=strand_regexmatch.group(1)


        if breakpoint_regexmatch:
            break

        result_row = list(result_dictionary.values())
        with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/"+ selection_result, "a", newline="") as add_result:
            add_result = csv.writer(add_result, delimiter="\t")
            add_result.writerow(result_row)



def obtain_scores():
    import csv
    with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/" + score_file, "w") as write_header:
        header = ['Mismatch Risk', 'Azimuth Score']
        write_header.write(header,delimiter="\t")


    with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/" + score_file, "a") as write_scores:
        with open("/home/zhang.r/CRISPR/dsnickfury/dsNickFury3PlusOrchid/selection_withelevation/" + selection_result, "r") as read_scores:
            for line in read_scores:
                write_scores.writerow(line[1],line[2])



gene_list = ["HPRT1","NF1", "NF2", "CCDC101", "MED12", "TADA2B","TADA1", "CD33", "CD13", "CD15","CUL3"]
chr_list=["X","17","22","16","X","4","1","19","15",'11',"2"]
start_list=["134460153","31094927","29603556","28553915","71118556","7041899","166856510","51225064","89784889","94543840","224470150"]
end_list=["134520513","31382116","29698598","28591790","71142454","7057952","166876327","51243860","89815401","94549898","224585397"]
for i in range (0,8):
    gene_name=gene_list[i]
    ensembl_chr=chr_list[i]
    ensembl_start=start_list[i]
    ensembl_end=end_list[i]
    fasta_file = gene_name + "_sequence.fa"
    selection_output = gene_name + "_output.txt"
    selection_result = gene_name + "_result.txt"
    score_file = gene_name + "_score.txt"
    get_genesequence()
    selection()
    result_file()
    selection__result()




