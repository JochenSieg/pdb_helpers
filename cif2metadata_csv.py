import os
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import gzip
from joblib import Parallel, delayed

def get_mmcif_dict(path: Path) -> dict:
    """
    Function that wraps biopythons MMCIF2Dict class for parsing data in mmCIF files to python dict.
    :param path: File path to the mmCIF file.
    :return: dict with CIF data if successful else None.
    """
    mmcif_path = Path(path)
    suffixes = mmcif_path.suffixes

    if len(suffixes) > 0:
        if suffixes[-1] == '.cif':
            with open(mmcif_path, 'r') as f:
                d = MMCIF2Dict(f)
        elif len(suffixes) > 1 and suffixes[-1] == '.gz' and suffixes[-2] == '.cif':
            with gzip.open(mmcif_path, 'rt') as f:
                d = MMCIF2Dict(f)
        else:
            return None
    else:
        print('ERROR: No cif format extension in {}'.format(mmcif_path))
        return None

    if '_entry.id' not in d.keys():
        print('ERROR: File seems not to be mmCIF format: {}'.format(mmcif_path))
        return None
    return d


def parse_seq(orginal_seq):
    # TAKEN FROM HHSUITE SOEDING ET AL. https://github.com/soedinglab/hh-suite
    # LICENSED UNDER GNU General Public License v3.0
    """Parses the cif fasta sequence and replaces non-canonical residues with their canonical counterparts."""

    THREE2ONE = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
        'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
        'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'MSE': 'M',
        'HYP': 'P', 'MLY': 'K', 'SEP': 'S', 'TPO': 'T', 'CSO': 'C', 'PTR': 'Y', 'KCX': 'K',
        'CME': 'C', 'CSD': 'A', 'CAS': 'C', 'MLE': 'L', 'DAL': 'A', 'CGU': 'E', 'DLE': 'L',
        'FME': 'M', 'DVA': 'V', 'OCS': 'C', 'DPR': 'P', 'MVA': 'V', 'TYS': 'Y', 'M3L': 'K',
        'SMC': 'C', 'ALY': 'K', 'CSX': 'C', 'DCY': 'C', 'NLE': 'L', 'DGL': 'E', 'DSN': 'S',
        'CSS': 'C', 'DLY': 'K', 'MLZ': 'K', 'DPN': 'F', 'DAR': 'R', 'PHI': 'F', 'IAS': 'D',
        'DAS': 'D', 'HIC': 'H', 'MP8': 'P', 'DTH': 'T', 'DIL': 'I', 'MEN': 'N', 'DTY': 'Y',
        'CXM': 'M', 'DGN': 'G', 'DTR': 'W', 'SAC': 'S', 'DSG': 'N', 'MME': 'M', 'MAA': 'A',
        'YOF': 'Y', 'FP9': 'P', 'FVA': 'V', 'MLU': 'L', 'OMY': 'Y', 'FGA': 'E', 'MEA': 'F',
        'CMH': 'C', 'DHI': 'H', 'SEC': 'C', 'OMZ': 'Y', 'SCY': 'C', 'MHO': 'M', 'MED': 'M',
        'CAF': 'C', 'NIY': 'Y', 'OAS': 'S', 'SCH': 'C', 'MK8': 'L', 'SME': 'M', 'LYZ': 'K'
    }

    CANONICAL_RESIDUES = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                          'Y'}

    seq = orginal_seq

    while seq.find('(') != -1:
        start_pos = seq.find('(')
        stop_pos = seq.find(')')
        residue = seq[start_pos + 1:stop_pos]

        try:
            canonical = THREE2ONE[residue]
        except KeyError:
            canonical = 'X'

        if start_pos == 0:
            seq = canonical + seq[stop_pos + 1:]
        elif stop_pos == len(seq):
            seq = seq[0:start_pos] + canonical
        else:
            pre_seq = seq[0:start_pos]
            post_seq = seq[stop_pos + 1:]
            seq = pre_seq + canonical + post_seq

    seq = seq.replace('\n', '')
    seq_array = []

    for c in seq:
        if c in CANONICAL_RESIDUES:
            seq_array.append(c)
        else:
            seq_array.append('X')

    seq = ''.join(seq_array)

    return seq


class SourceOrganismHandler:
    def __init__(self, source_organism_data: dict):
        """
        Handles different aspects of the source organism annotation of a single PDB entry.
        :param source_organism_data: dict of source organism fields in the CIF format
        """
        self.data_dict = source_organism_data

        # We map PDBid to an organism list.
        #   -If pdbx_organism_scientific is defined we use this name else pdbx_gene_src_scientific_name
        #   (following soeding et al. in https://github.com/soedinglab/hh-suite/blob/master/scripts/cif2fasta.py).
        #   -If neither annotation is provided we drop the PDBid completely.
        if self.data_dict['_entity_src_nat.pdbx_organism_scientific'] is not None:
            self.organism_list = self.data_dict['_entity_src_nat.pdbx_organism_scientific']
        elif self.data_dict['_entity_src_gen.pdbx_gene_src_scientific_name'] is not None:
            self.organism_list = self.data_dict['_entity_src_gen.pdbx_gene_src_scientific_name']
        else:
            self.organism_list = []

    def is_chimeric(self):
        """
        Checks if all entities in the PDB/CIF entry are from the same source organism.
        NOTE: Obviously this function is just as good as the annotation itself. Typos and other annotation errors
              will lead to false positive results.
        :return: True if different source organisms are in the entry,
                 False if only a single source organism is in the entry,
                 None if no source organism is in the entry.
        """
        organism_set = set([org.strip().lower() for org in self.organism_list])

        # empty lists == no species annotated -> ignore
        if len(organism_set) < 1:
            return None

        # remove entries which are just subsets of other entries
        # this is often the case, for example there is the organism with strain and without strain.
        org_list = list(organism_set)
        organism_set = set()
        for org1 in org_list:
            is_subset = False
            for org2 in org_list:
                if org1 != org2:
                    if org1 in org2:
                        is_subset = True
                        break
            if not is_subset:
                organism_set.add(org1)

        if len(organism_set) > 1:
            # is chimeric
            return True

        return False


class CIFMetaDataExtractor:
    def __init__(self, path: Path, min_seq_len=30, skip_non_standard_only=False):
        """
        Creates an extractor instance.
        :param path: Path to CIF file.
        :param min_seq_len: Minimal length of a protein sequence to consider it.
        :param skip_non_standard_only: Whether to ignore protein chains consisting only of non-standard amino acids.
        """
        self.min_seq_len = min_seq_len
        self.skip_non_standard_only = skip_non_standard_only
        self.d = get_mmcif_dict(path)
        self.pdbid = self.d['_entry.id']
        if len(self.pdbid) != 1:
            print('WARNING: No or multiple PDB IDs: {} in file: {}'.format(self.pdbid, path))
        self.pdbid = self.pdbid[0]

    def is_valid(self):
        """
        Returns if CIF file/entry is valid.
        :return: True if valid else false.
        """
        return self.d is not None

    @staticmethod
    def get_header_fields():
        """
        Returns the header describing the metadata fields returned by get_meta_data_items.
        :return: list of field names.
        """
        fields = ['pdbchainid', 'seq_len', 'is_chimeric_entity', 'resolution', 'rfree', 'exp_method', 'src_org',
                  'is_chimeric_entry', 'contains_dna', 'contains_rna', 'ec_str']
        return fields

    def get_meta_data_items(self):
        """
        Main function to retrieve metadata annotation for all protein chain entities in the CIF entry.
        :return: yields meta data info for each chain as tuple.
        """

        if '_entity_poly.entity_id' not in self.d:
            print('ERROR: {} has no poly entity entry: _entity_poly.entity_id'.format(self.pdbid))
            return

        # First, we get global metadata of the whole PDB entry

        resolution = self.get_resolution()
        rfree = self.get_rfree()
        exp_method = self.get_experimental_method()

        source_organisms_data = self.get_source_organims()
        is_chimeric_entry = SourceOrganismHandler(source_organisms_data).is_chimeric()

        contains_dna = 'polydeoxyribonucleotide' in self.d['_entity_poly.type']
        contains_rna = 'polyribonucleotide' in self.d['_entity_poly.type']

        # Second, we iterate over each polypeptide chain
        peptide_entity_ids = self.d['_entity_poly.entity_id']

        # iterate over entity_poly
        for i, poly_id in enumerate(peptide_entity_ids):

            # polypeptide chain
            if self.d['_entity_poly.type'][i] == 'polypeptide(L)':
                chainids = self.d['_entity_poly.pdbx_strand_id'][i]
                seq = parse_seq(self.d['_entity_poly.pdbx_seq_one_letter_code'][i])

                if len(seq) < self.min_seq_len:
                    # skip too short chains
                    continue

                if self.skip_non_standard_only:
                    # remove chains that contain only unknown residues
                    reduced_set = set(seq)
                    if len(reduced_set) == 1 and next(iter(reduced_set)) == 'X':
                        continue

                src_org, is_chimeric_entity = self.get_chain_source_organisms(poly_id)
                ec_str = self.get_chain_ec_numbers(poly_id)

                # iterate over identical chains (e.g. homo-mers)
                for chain in chainids.split(','):
                    pdbchainid = '{}_{}'.format(self.pdbid, chain)
                    yield pdbchainid, len(seq), is_chimeric_entity,\
                          resolution, rfree, exp_method, src_org,\
                          is_chimeric_entry, contains_dna, contains_rna, ec_str,
            else:
                # non peptide chain
                pass

    def get_chain_ec_numbers(self, poly_id):
        """
        Get EC-numbers for this chain.
        :param poly_id: The poly_id for this entity.
        :return: Comma separated string of all unique EC-numbers. None if no EC-Number available.
        """
        ec_set = set()
        if '_entity.pdbx_ec' in self.d.keys():
            for i, ecnum in enumerate(self.d['_entity.pdbx_ec']):
                if self.d['_entity.id'][i] == poly_id:
                    ec_set.add(ecnum)
        elif '_entity_keywords.ndb_ec' in self.d.keys():
            for i, ecnum in enumerate(self.d['_entity_keywords.ndb_ec']):
                if self.d['_entity.id'][i] == poly_id:
                    ec_set.add(ecnum)

        if len(ec_set) == 1:
            # we got one EC number. If it is a missing value as indicated by '?' we return None.
            ec_str = next(iter(ec_set))
            return None if ec_str == '?' else ec_str
        elif len(ec_set) < 1:
            # no EC number. Return None.
            return None
        else:
            # multiple EC numbers for this chain. Create comma separated list from all
            # TODO what about '?'
            return ','.join(sorted(ec_set))

    def get_chain_source_organisms(self, poly_id: str):
        """
        Get source organisms of the chain corresponding to the poly_id.
        :param poly_id: The identifier for this entity.
        :return: pair of a comma separated string of all source organisms and a bool indicating whether the chain
                 is chimeric (is from multiple source organims).
        """
        is_chimeric_entity = False
        chain_src_org_set = set()

        # Select the source organism for the chain entity. Prefer _entity_src_nat entry.
        if '_entity_src_nat.pdbx_organism_scientific' in self.d.keys():
            for j, entity_id in enumerate(self.d['_entity_src_nat.entity_id']):
                if entity_id == poly_id:
                    chain_src_org_set.add(self.d['_entity_src_nat.pdbx_organism_scientific'][j].strip())
        elif '_entity_src_gen.pdbx_gene_src_scientific_name' in self.d.keys():
            for j, entity_id in enumerate(self.d['_entity_src_gen.entity_id']):
                if entity_id == poly_id:
                    chain_src_org_set.add(self.d['_entity_src_gen.pdbx_gene_src_scientific_name'][j].strip())

        # Differentiate cases: No source organism, multiple source organisms (chimeric chain) and one-2-one.
        if len(chain_src_org_set) < 1:
            # This chain has no organism annotated
            is_chimeric_entity = False
            src_org = None
        elif len(chain_src_org_set) > 1:
            # this chain is the chimeric product from multiple organisms.
            # We annotated that with a comma separated list of the organisms.
            src_org = ','.join(sorted(chain_src_org_set))
            is_chimeric_entity = True
        else:
            # The chain is from exactly one organism
            is_chimeric_entity = False
            src_org = next(iter(chain_src_org_set))

        return src_org, is_chimeric_entity

    def get_resolution(self):
        """
        Get resolution of the cif entry.
        :return: float or None
        """
        resolution_fields = ['_refine.ls_d_res_high', '_reflns.d_resolution_high', '_em_3d_reconstruction.resolution']

        for field in resolution_fields:
            if field in self.d.keys():
                reso_list = self.d[field]
                if len(reso_list) != 1:
                    print('ERROR: {} Unexpected number of resolution: {}'.format(self.pdbid, reso_list))
                try:
                    return float(reso_list[0])
                except ValueError:
                    return None
        return None

    def get_experimental_method(self):
        """
        Get experimental method if available.
        :return: Experimental method if availabel or None.
        """
        if '_exptl.method' in self.d.keys():
            method_list = self.d['_exptl.method']
            if len(method_list) != 1:
                # TODO one might argue if multiple methods are an error.
                print('ERROR: {} Unexpected number of experimental method: {}'.format(self.pdbid, method_list))
            return method_list[0].strip().replace('\n', ' ')
        return None

    def get_rfree(self):
        """
        Get rfree value if available.
        :return: float if present else None
        """
        if '_refine.ls_R_factor_R_free' in self.d.keys():
            rfree_list = self.d['_refine.ls_R_factor_R_free']
            if len(rfree_list) != 1:
                print('ERROR: {} Unexpected number of rfree: {}'.format(self.pdbid, rfree_list))
            try:
                return float(rfree_list[0])
            except ValueError:
                return None
        return None

    def get_source_organims(self):
        """
        Returns source organisms info of the mmcif entry
        :return: dict with relevant fields.
        """
        org_sci = None
        gene_src = None
        ncbi_id = None

        # Required in PDB entries: no
        # Required for PDB deposition: yes
        # Used in current PDB entries: Yes, in about 10.1 % of entries
        if '_entity_src_nat.pdbx_organism_scientific' in self.d.keys():
            org_sci = self.d['_entity_src_nat.pdbx_organism_scientific']

        # Required in PDB entries: no
        # Required for PDB deposition: yes
        # Used in current PDB entries: Yes, in about 87.6 % of entries
        if '_entity_src_gen.pdbx_gene_src_scientific_name' in self.d.keys():
            gene_src = self.d['_entity_src_gen.pdbx_gene_src_scientific_name']

        # Required in PDB entries: no
        # Required for PDB deposition: yes
        # Used in current PDB entries: Yes, in about 10.0 % of entries
        if '_entity_src_nat.pdbx_ncbi_taxonomy_id' in self.d.keys():
            ncbi_id = self.d['_entity_src_nat.pdbx_ncbi_taxonomy_id']

        return {'_entity_src_nat.pdbx_organism_scientific': org_sci,
                '_entity_src_gen.pdbx_gene_src_scientific_name': gene_src,
                '_entity_src_nat.pdbx_ncbi_taxonomy_id': ncbi_id
                }


def get_all_cif_files(directory: Path):
    """
    Returns all *.cif/*.cif.gz paths in a directory.
    :param directory: The directory.
    :return: List of Paths to CIF files.
    """
    allowed_file_extensions = ['.cif', '.cif.gz']
    return [Path(entry.path) for entry in os.scandir(directory)
            if any(entry.name.endswith(ext) for ext in allowed_file_extensions)]


def run(cif_path: Path):
    """
    Wrapper function to execute CIF metadata parsing in parallel.
    :param cif_path: Path to the CIF file.
    :return: list of metadata. One entry in the list describes the meta data of a protein chain entity.
    """
    extractor = CIFMetaDataExtractor(cif_path)
    if extractor.is_valid():
        return list(extractor.get_meta_data_items())
    return list()


def write_tsv(header, entry_rows: list, outfile: Path, additional_info: dict):
    """
    Writes a tsv file with the metadata to disc. One row corresponds to one protein chain entity.
    :param header:
    :param entry_rows:
    :param outfile:
    :param additional_info:
    :return:
    """
    sep = '\t'
    add_info_keys = sorted(additional_info.keys())
    i = 0
    with open(outfile, 'w') as f:

        # write header
        f.write(sep.join(header) + '\n')

        # loop over PDB/CIF entries
        for rows in entry_rows:
            # loop over chains
            for row in rows:
                f.write(sep.join(map(str, row)))

                for k in add_info_keys:
                    f.write('{}{}'.format(sep, additional_info[k][i]))

                f.write('\n')
                i += 1


def read_pfam_pdbmapping(path: Path) -> dict:
    """
    Read the Pfam PDB chain mapping file from the Pfam ftp server into memory.
    :param path: Path to the mapping file.
    :return: A look up dict that maps PDB chain IDs to their PfamIDs.
    """
    def parse_pfam(filehandle):
        pfam_dict = {}  # maps PDBchainid to PfamID

        for line in filehandle:
            cols = line.split(';')
            pdbid = cols[0].strip()
            chainid = cols[1].strip()
            pfamid = cols[3].strip()

            identifier = '{}_{}'.format(pdbid, chainid)

            if identifier not in pfam_dict.keys():
                pfam_dict[identifier] = []
            pfam_dict[identifier].append(pfamid)

        return pfam_dict

    if path.name.endswith('.gz'):
        with gzip.open(path, 'rt') as f:
            pfam_dict = parse_pfam(f)
    else:
        with open(path, 'r') as f:
            pfam_dict = parse_pfam(f)

    return pfam_dict


def get_pfam_id(pfam_dict, entry_rows):
    """
    Collects Pfam identifiers for PDBchainIDs in entry_rows.
    :param pfam_dict: The lookup dict for PDBchainID to PfamID.
    :param entry_rows: The data to make look ups for.
    :return: A list of PfamIDs (corresponding to the flattened list of PDBchainIDs in entry_rows).
    """
    pfam_ids = []
    if pfam_dict:
        # loop over PDB/CIF entries
        for rows in entry_rows:
            # loop over chains
            for row in rows:
                pdbchainid = row[0]
                if pdbchainid in pfam_dict.keys():
                    # There are pfam ids for this chain
                    pfam_ids.append(','.join(pfam_dict[pdbchainid]))
                else:
                    pfam_ids.append(None)
    return pfam_ids


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description=
        """
        Reads all *.cif/*.cif.gz files in a directory, extracts metadata information and writes this information to
        a CSV. Metadata is extracted on the basis of polypeptide chains. This means there is one row in the resulting
        CSV for each protein chain. DNA/RNA and other non-protein-chains will be ignored.
        
        The metadata extracted will describe information on the protein chain entity, for example cahin length. In 
        addition information of the whole PDB entry is included, for example the resolution, rfree or if the entry
        contains any DNA.
        
        More information on the schema for molecular entities in mmCIF of the PDB can be found here:
        https://mmcif.wwpdb.org/docs/tutorials/content/molecular-entities.html
        
        
        Example:
        
        python cif2metadata_csv.py --directory /data/pdb/current/data/structures/all/mmCIF/ --outfile pdb_metadata.csv
               --pfam /local/sieg/pfam/pdbmap.gz -j 6
        
        """)

    parser.add_argument('--directory', '-d', required=True, type=str,
                        help='Directory of *.cif files. Can also be *.cif.gz.')
    parser.add_argument('--outfile', '-o', required=True, type=str,
                        help='Output CSV file.')
    parser.add_argument('--jobs', '-j', type=int, default=1,
                        help='Number of jobs for parallel execution.')
    parser.add_argument('--pfam', type=str, default=None,
                        help='Path to Pfam PDB mapping file. File is called pdbmap on Pfam ftp server.')

    args = parser.parse_args()

    input_dir = Path(args.directory)
    outfile = Path(args.outfile)
    jobs = args.jobs
    pfam_file = args.pfam

    if not input_dir.is_dir():
        print('Error: CSV file does not exist.')
        sys.exit(1)

    # get all cif files in directory
    cifs = get_all_cif_files(input_dir)

    pfam_dict = read_pfam_pdbmapping(Path(pfam_file)) if pfam_file else None

    with Parallel(n_jobs=jobs) as parallel:
        entry_rows = parallel(delayed(run)(c) for c in cifs)

    header = CIFMetaDataExtractor.get_header_fields()

    # add additional info
    additional_info = {}
    if pfam_dict:
        additional_info['pfam_id'] = get_pfam_id(pfam_dict, entry_rows)
        header.append('pfam_id')

    write_tsv(header, entry_rows, outfile, additional_info)


if __name__ == "__main__":
    main()
