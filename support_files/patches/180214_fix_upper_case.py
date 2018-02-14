# run within the Qiita environment, suggest using a screen session

from qiita_db.study import Study
from qiita_db.software import Software
from biom import load_table
from qiita_db.util import compute_checksum
from qiita_db.sql_connection import TRN
from biom.util import biom_open
from os import rename
from collections import defaultdict


# this is necessary to overcome this issue
# https://github.com/biocore/biom-format/issues/761
def collapse_f(table, axis):
    return table.sum(axis=axis)


studies = Study.get_by_status('private').union(
    Study.get_by_status('public')).union(Study.get_by_status('sandbox'))

sft = Software.from_name_and_version('deblur', '1.0.3')
# [0] deblur only has 1 command
cmd = sft.commands[0]

artifacts = [a for s in studies for a in s.artifacts()
             if a.processing_parameters is not None
             and a.processing_parameters.command == cmd]


sql = "UPDATE qiita.filepath SET checksum = %s WHERE filepath_id = %s"
for a in artifacts:
    # putting all this in a transaction in case something fails it does
    # it nicely
    with TRN:
        ftps = {ft: (fid, fp) for fid, fp, ft in a.filepaths
                if ft in ['biom', 'preprocessed_fasta']}
        biom = ftps['biom'][1]
        t = load_table(biom)
        current = set(t.ids('observation'))
        updated = set(map(lambda x: x.upper(), current))
        difference = current ^ updated
        if difference:
            print ('*********>\nChanging biom: ', a.id, biom)
            # checking for duplicated ids
            if len(current) != len(updated):
                duplicates = defaultdict(list)
                # getting the main list
                for key in difference:
                    if key in current:
                        duplicates[key.upper()].append(key)
                # adding cases where the key is in the biom in an upper form
                for key in duplicates.keys():
                    if key in current:
                        duplicates[key].append(key)

                # formatting for easier processing
                to_merge = {}
                ids_to_replace = {}
                for k, v in duplicates.items():
                    # 1 means that is the regular upper change
                    if len(v) == 1:
                        ids_to_replace[v[0]] = k
                    else:
                        for vv in v:
                            to_merge[vv] = k
                merge_fn = (lambda id_, x: to_merge[id_]
                            if id_ in to_merge else id_)
                t = t.collapse(merge_fn, norm=False, min_group_size=1,
                               axis='observation', collapse_f=collapse_f)
            else:
                ids_to_replace = {c: c.upper() for c in current
                                  if c != c.upper()}

            t.update_ids(ids_to_replace, axis='observation', strict=False,
                         inplace=True)

            with biom_open(biom, 'w') as f:
                t.to_hdf5(f, t.generated_by)
            checksum = compute_checksum(biom)

            TRN.add(sql, [checksum, ftps['biom'][0]])

            fna = ftps['preprocessed_fasta'][1]
            tmp = fna + '.tmp'
            with open(tmp, 'w') as out:
                for seq in t.ids('observation'):
                    out.write('>%s\n%s\n' % (seq, seq))
            rename(tmp, fna)
            checksum = compute_checksum(fna)

            TRN.add(sql, [checksum, ftps['preprocessed_fasta'][0]])

            TRN.execute()
