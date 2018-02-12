# run within the Qiita environment, suggest using a screen session

from qiita_db.study import Study
from qiita_db.software import Software
from biom import load_table
from qiita_db.util import compute_checksum
from qiita_db.sql_connection import TRN
from biom.util import biom_open
from skbio.io import read
from os import rename, remove


studies = Study.get_by_status('private').union(
    Study.get_by_status('public')).union(Study.get_by_status('sandbox'))

sft = Software.from_name_and_version('deblur', '1.0.4')
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
        current = t.ids('observation')
        updated = map(lambda x: x.upper(), current)
        if set(current) ^ set(updated):
            print '*********>\nChanging biom: ', a.id, fp
                t.update_ids({i: i.upper() for i in t.ids('observation')},
                             axis='observation', inplace=True)
            with biom_open(fp, 'w') as f:
                t.to_hdf5(f, t.generated_by)
            checksum = compute_checksum(fp)

            TRN.add(sql, [checksum, ftps['biom'][0]])

            fna = ftps['preprocessed_fasta'][1]
            tmp = fna + '.tmp'
            with open(tmp, 'w') as out:
                for seq in read(fna, format='fasta'):
                    seq = str(seq)
                    sequ = seq.upper()
                    out.write('>%s\n%s\n' % (sequ, sequ))
            rename(tmp, fna)
            checksum = compute_checksum(fna)

            TRN.add(sql, [checksum, ftps['preprocessed_fasta'][0]])
