#!/usr/bin/python3

# API Reference https://ncbi.github.io/blast-cloud/dev/api.html
# code copied, adapted and simplified from https://github.com/kpodkalicki/BLAST-API-Implementation

import re
import time
import shutil
import requests
from datetime import datetime

_API_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

_Delta_Wait = 5
_CBSs = [0, 1, 2, 3]
_NCBI_GIs = ['T', 'F']
_Filter = ['F', 'T', 'L', 'mT', 'mL']
_Programs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
_Format_Types = ['HTML', 'Text', 'XML', 'XML2', 'JSON2', 'Tabular']
_Scoring_Matrices = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM250', 'PAM30', 'PAM70']


def search(query: str, database: str, program: str, additional_params: dict = None,
           expect: int = None,  threshold: int = None, word_size: int = None,
           alignments: int = None, num_threads: int = None, nucl_reward: int = None,
           nucl_penalty: int = None, hitlist_size: int = None, descriptions: int = None,
           filter: str = None, ncbi_gi: str = None, matrix: str = None, format_type: str = None, cbs: str = None,
           gapcosts: (int, int) = None) -> (str, int):
    r"""Sends search submission to NCBI-BLAST Common URL API.

    :param query: Search query. Accession, GI, or FASTA.
    :param database: Name of existing database or one uploaded to blastdb_custom.
    :param program: BLAST Program. One of: ['blastn', 'megablast', 'blastp', 'blastx', 'tblastn', 'tblastx'].
    :param filter: Low complexity filtering. F to disable. T or L to enable. Prepend “m” for mask at lookup (e.g., mL).
    :param format_type: Report type. One of: ['HTML', 'Text', 'XML', 'XML2', 'JSON2', 'Tabular']. Default: 'HTML'.
    :param expect: Expect value. Number greater than zero.
    :param nucl_reward: Reward for matching bases (BLASTN and megaBLAST). Integer greater than zero.
    :param gapcosts: Gap existence and extension costs. Tuple of two positive integers.
    :param matrix: Scoring matrix name. Default: 'BLOSUM62' One of:
            ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM250', 'PAM30' or 'PAM70'].
    :param hitlist_size: Number of databases sequences to keep. Integer greater than zero.
    :param descriptions: Number of descriptions to print (applies to HTML and Text). Integer greater than zero.
    :param alignments: Number of alignments to print (applies to HTML and Text). Integer greater than zero.
    :param ncbi_gi: Show NCBI GIs in report. 'T' or 'F'
    :param threshold: Neighboring score for initial words. Positive integer (BLASTP default is 11).
            Does not apply to BLASTN or MegaBLAST.
    :param word_size: Size of word for initial matches. Positive integer.
    :param cbs: Composition based statistics algorithm to use. One of [0, 1, 2, 3].
            See comp_based_stats in https://www.ncbi.nlm.nih.gov/books/NBK279684/ for details.
    :param num_threads: Number of virtual CPUs to use. 	Integer greater than zero (default is 1).
    :return: Tuple of request_id and estimated_time in seconds until the search is completed.
    """

    if program not in _Programs:
        print(f"program {program} not supported")
        print(f"Available programs: {_Programs}")
        exit(1)

    params = {
        'CMD': 'Put',
        'QUERY': query,
        'DATABASE': database,
        'PROGRAM': program
    }

    # combining with additional parameters
    if additional_params is not None:
        params = dict(list(additional_params.items())+list(params.items()))

    # Optional Params

    # FORMAT_OBJECT	Object type	String	Get	SearchInfo (status check) or Alignment (report formatting).

    if expect is not None:
        if expect > 0: params['EXPECT'] = str(expect)
        else: print(f"expect {expect} not greater than 0"); exit(2)
    if threshold is not None:
        if threshold >= 0: params['THRESHOLD'] = threshold
        else: print(f"threshold {threshold} not greater than 0"); exit(2)
    if word_size is not None:
        if word_size >= 0: params['WORD_SIZE'] = word_size
        else: print(f"word_size {word_size} not greater than 0"); exit(2)
    if alignments is not None:
        if alignments > 0: params['ALIGNMENTS'] = str(alignments)
        else: print(f"alignments {alignments} not greater than 0"); exit(2)
    if num_threads is not None:
        if num_threads > 0: params['NUM_THREADS'] = num_threads
        else: print(f"num_threads {num_threads} not greater than 0"); exit(2)
    if nucl_reward is not None:
        if nucl_reward > 0: params['NUCL_REWARD'] = str(nucl_reward)
        else: print(f"nucl_reward {nucl_reward} not greater than 0"); exit(2)
    if nucl_penalty is not None:
        if nucl_reward > 0: params['NUCL_PENALTY'] = str(nucl_penalty)
        else: print(f"nucl_penalty {nucl_penalty} not greater than 0"); exit(2)
    if hitlist_size is not None:
        if hitlist_size > 0: params['HITLIST_SIZE'] = str(hitlist_size)
        else: print(f"hitlist_size {hitlist_size} not greater than 0"); exit(2)
    if descriptions is not None:
        if descriptions > 0: params['DESCRIPTIONS'] = str(descriptions)
        else: print(f"descriptions {descriptions} not greater than 0"); exit(2)

    if filter is not None:
        if filter in _Filter: params['FILTER'] = filter
        else: print(f"filter {filter} not in {_Filter}"); exit(3)
    if ncbi_gi is not None:
        if ncbi_gi in _NCBI_GIs: params['NCBI_GI'] = ncbi_gi
        else: print(f"filter {ncbi_gi} not in {_NCBI_GIs}"); exit(3)
    if matrix is not None:
        if matrix in _Scoring_Matrices: params['MATRIX'] = matrix
        else: print(f"filter {matrix} not in {_Scoring_Matrices}"); exit(3)
    if format_type is not None:
        if format_type in _Format_Types: params['FORMAT_TYPE'] = format_type
        else: print(f"filter {format_type} not in {_Format_Types}"); exit(3)
    if cbs is not None:
        if cbs in _CBSs: params['COMPOSITION_BASED_STATISTICS'] = str(cbs)
        else: print(f"filter {cbs} not in {_CBSs}"); exit(3)

    if gapcosts is not None:
        existence, extension = gapcosts
        if existence > 0 and extension > 0: params['GAPCOSTS'] = f"{existence} {extension}"
        else: print(f"gapcosts {gapcosts}: either one or both is geater then 0"); exit(2)

    print(f"Params: {params}")

    response = requests.get(_API_URL, params)

    # print(response.text)

    qblast_info = __cropp_qblast_info__(response.text)
    return qblast_info['RID'], int(qblast_info['RTOE'])


def get_results(request_id, format_type='HTML', hitlist_size=None, descriptions=None, alignments=None,
                ncbi_gi=None, format_object=None, results_file_path='results.zip'):
    r"""Retrieves results from NCBI.

    :param request_id: ID of requested submission.
    :param format_type: Report type. One of: ['HTML', 'Text', 'XML', 'XML2', 'JSON2', 'Tabular']. Default: 'HTML'.
    :param hitlist_size: Number of databases sequences to keep. Integer greater than zero.
    :param descriptions: Number of descriptions to print (applies to HTML and Text). Integer greater than zero.
    :param alignments: Number of alignments to print (applies to HTML and Text). Integer greater than zero.
    :param ncbi_gi: Show NCBI GIs in report. 'T' or 'F'
    :param format_object: Object type. SearchInfo (status check) or Alignment (report formatting).
            Only Alignment is valid for retrieving results.
    :param results_file_path: Results relative file path (applies to XML2 and JSON2).
    :return: Results or relative path to results file.
    """

    params = {
        'CMD': 'Get',
        'RID': request_id,
        'FORMAT': format_type,
    }

    if hitlist_size is not None:  params['HITLIST_SIZE'] = hitlist_size
    if descriptions is not None:  params['DESCRIPTIONS'] = descriptions
    if alignments is not None:    params['ALIGNMENTS'] = alignments
    if ncbi_gi is not None:       params['NCBI_GI'] = ncbi_gi
    if format_object is not None: params['FORMAT_OBJECT'] = format_object

    if format_type == 'XML2' or format_type == 'JSON2':
        response = requests.get(_API_URL, params, stream=True)
        with open(results_file_path, 'wb') as file:
            shutil.copyfileobj(response.raw, file)
        return results_file_path
    else:
        response = requests.get(_API_URL, params)
        return response.text


def wait_for_results(request_id, estimated_time=2):
    r"""Waits for availability of results and retrieves it.

    :param request_id: ID of requested submission
    :param estimated_time: estimated time in seconds until the search is completed
    :return: Results or relative path to results file
    """

    time.sleep(int(estimated_time))
    status = __check_submission_status(request_id)
    while status == 'WAITING':
        time.sleep(_Delta_Wait)
        status = __check_submission_status(request_id)
        print(f"Submission status: {status} {datetime.now().strftime('%H:%M:%S')}")

    if status == 'UNKNOWN':
        raise ValueError("NCBI returned 'UNKNOWN' submission status")

    return status


def __check_submission_status(request_id):
    r"""Checks submission status.

    :param request_id: ID of requested submission
    :return: Status of submission. 'WAITING', 'UNKNOWN' or 'READY'
    """

    params = {
        'CMD': 'Get',
        'FORMAT_OBJECT': 'SearchInfo',
        'RID': request_id
    }

    response = requests.get(_API_URL, params)
    return __cropp_qblast_info__(response.text)['Status']


def __cropp_qblast_info__(html):  # why the f*ck did he parse html with regex when there are other output formats?!
    search_results = re.findall('QBlastInfoBegin(.*?)QBlastInfoEnd', html, flags=re.DOTALL)
    if search_results:
        result_dict = dict()
        for search_result in search_results:
            search_result = search_result.strip()
            entries = search_result.split("\n")
            for entry in entries:
                key_value_pair = entry.split('=')
                if len(key_value_pair) > 0:
                    key = key_value_pair[0].strip()
                    value = None
                    if len(key_value_pair) > 1:
                        value = key_value_pair[1].strip()
                    result_dict[key] = value
        return result_dict
    return dict()


if __name__ == '__main__':
    request_id, estimated_time = search(query="NM_001260153.1", database="nt", program="blastn",
                                        additional_params={'CLIENT': 'web'})
    print(f"Request id: {request_id} Estimated Time: {estimated_time}")

    print(wait_for_results(request_id, estimated_time))

    print('#########################  HTML  #########################')
    print(get_results(request_id, format_type='HTML'))
    print('\n\n\n')

    print('#########################  Text  #########################')
    print(get_results(request_id, format_type='Text'))
    print('\n\n\n')

    print('#########################  XML  #########################')
    print(get_results(request_id, format_type='XML'))
    print('\n\n\n')

    print(get_results(request_id, format_type='XML2', results_file_path='ex1_xml2.zip'))
    print(get_results(request_id, format_type='JSON2', results_file_path='ex1_json2.zip'))

    print('#########################  XML  #########################')
    print(get_results(request_id, format_type='Tabular'))
    print('\n\n\n')
