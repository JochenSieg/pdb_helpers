import json
import requests
import pandas as pd

"""
This file shows how to communicate with the PDBs GraphQL API via Python.

Use the get_meta_info_for_pdbids function to get a DataFrame of metadata information for a list of PDB ids.
"""


PDB_GRAPHQL = 'https://data.rcsb.org/graphql'


def get_meta_info_query(pdbids: list) -> dict:
    return {
        'query': """{{
          entries(entry_ids: {}) {{
            rcsb_id
            exptl {{
              method
            }}
            rcsb_entry_info {{
              polymer_entity_count_DNA
              polymer_entity_count_RNA
              polymer_entity_count_nucleic_acid
              polymer_entity_count_nucleic_acid_hybrid
              polymer_entity_taxonomy_count
              experimental_method_count
            }}
            pdbx_vrpt_summary {{
              PDB_resolution
            }}
          }}
        }}""".format(json.dumps(pdbids))
    }


def do_request(url: str, query: dict) -> dict:
    response = requests.post(url=url, json=query)

    # print(response.text)

    print('Statuscode:', response.status_code)

    if response.status_code == 200:
        return json.loads(response.content.decode('utf-8'))
    else:
        return {}


def parse_meta_info(data_dict: dict, result_dict: dict):
    na = pd.NA

    if 'rcsb_id' not in data_dict:
        raise ValueError('ERROR: No PDB. Meta info response corrupt.')

    result_dict['pdbid'].append(data_dict['rcsb_id'])
    #     print(data_dict['rcsb_id'])

    if 'exptl' in data_dict:
        result_dict['exptl_method'].append(
            [m_dict['method'] for m_dict in data_dict['exptl'] if 'method' in m_dict.keys()])
    else:
        result_dict['exptl_method'].append(na)

    if 'rcsb_entry_info' in data_dict:

        if 'experimental_method_count' in data_dict['rcsb_entry_info']:
            result_dict['exptl_method_count'].append(data_dict['rcsb_entry_info']['experimental_method_count'])
        else:
            result_dict['exptl_method_count'].append(na)

        if 'polymer_entity_count_DNA' in data_dict['rcsb_entry_info']:
            result_dict['polymer_entity_count_DNA'].append(data_dict['rcsb_entry_info']['polymer_entity_count_DNA'])
        else:
            result_dict['polymer_entity_count_DNA'].append(na)

        if 'polymer_entity_count_RNA' in data_dict['rcsb_entry_info']:
            result_dict['polymer_entity_count_RNA'].append(data_dict['rcsb_entry_info']['polymer_entity_count_RNA'])
        else:
            result_dict['polymer_entity_count_RNA'].append(na)

        if 'polymer_entity_count_nucleic_acid' in data_dict['rcsb_entry_info']:
            result_dict['polymer_entity_count_nucleic_acid'].append(
                data_dict['rcsb_entry_info']['polymer_entity_count_nucleic_acid'])
        else:
            result_dict['polymer_entity_count_nucleic_acid'].append(na)

        if 'polymer_entity_count_nucleic_acid_hybrid' in data_dict['rcsb_entry_info']:
            result_dict['polymer_entity_count_nucleic_acid_hybrid'].append(
                data_dict['rcsb_entry_info']['polymer_entity_count_nucleic_acid_hybrid'])
        else:
            result_dict['polymer_entity_count_nucleic_acid_hybrid'].append(na)

        if 'polymer_entity_taxonomy_count' in data_dict['rcsb_entry_info']:
            result_dict['polymer_entity_taxonomy_count'].append(
                data_dict['rcsb_entry_info']['polymer_entity_taxonomy_count'])
        else:
            result_dict['polymer_entity_taxonomy_count'].append(na)
    else:
        result_dict['exptl_method_count'].append(na)
        result_dict['polymer_entity_count_DNA'].append(na)
        result_dict['polymer_entity_count_RNA'].append(na)
        result_dict['polymer_entity_count_nucleic_acid'].append(na)
        result_dict['polymer_entity_count_nucleic_acid_hybrid'].append(na)
        result_dict['polymer_entity_taxonomy_count'].append(na)

    if 'pdbx_vrpt_summary' in data_dict:

        if data_dict['pdbx_vrpt_summary'] is not None:
            if 'PDB_resolution' in data_dict['pdbx_vrpt_summary']:
                result_dict['resolution'].append(data_dict['pdbx_vrpt_summary']['PDB_resolution'])
            else:
                result_dict['resolution'].append(na)
        else:
            # i do not know why but this strange case appears, for example for 1ZY8
            result_dict['resolution'].append(na)
    else:
        result_dict['resolution'].append(na)

    if len(result_dict['exptl_method'][-1]) > 1 or result_dict['exptl_method_count'][-1] > 1:
        print('DEBUG INFO:', result_dict['pdbid'][-1], result_dict['exptl_method'][-1],
              result_dict['exptl_method_count'][-1])


def parse_json_response(json_response: dict, parse_f) -> pd.DataFrame:
    if 'data' not in json_response.keys():
        raise ValueError('ERROR: Invalid json_response. No data found.')
    if 'entries' not in json_response['data'].keys():
        raise ValueError('ERROR: Invalid json_response. No entries in data found.')

    # falls irgendwann notwendig kann der meta_info kram als eigene Klasse abstrahiert werden
    fields = ['pdbid', 'exptl_method', 'exptl_method_count', 'polymer_entity_count_DNA',
              'polymer_entity_count_RNA', 'polymer_entity_count_nucleic_acid',
              'polymer_entity_count_nucleic_acid_hybrid', 'polymer_entity_taxonomy_count', 'resolution']
    result_dict = {f: [] for f in fields}

    for data_dict in json_response['data']['entries']:
        parse_f(data_dict, result_dict)

    return pd.DataFrame(result_dict)


def get_meta_info_for_pdbids(pdbids: list) -> pd.DataFrame:
    chunk_size = 10000
    dfs = []
    for i in range(0, len(pdbids), chunk_size):
        chunk = pdbids[i:i + chunk_size]
        query = get_meta_info_query(chunk)
        json_response = do_request(PDB_GRAPHQL, query)

        df = parse_json_response(json_response, parse_meta_info)
        dfs.append(df)

    df = pd.concat(dfs)
    return df
