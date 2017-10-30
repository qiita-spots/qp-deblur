-- moving deblur v1.0.2 artifacts to v1.0.3
-- the parameters don't need to be changed as they are exactly the same

DO $do$
DECLARE
    old_deblur_id BIGINT;
    new_deblur_id BIGINT;
    a_info        RECORD;
    pj_id         UUID;
BEGIN
    SELECT command_id INTO old_deblur_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'deblur' AND s.version = '1.0.2' AND sc.name = 'deblur-workflow';

    SELECT command_id INTO new_deblur_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'deblur' AND s.version = '1.0.3' AND sc.name = 'deblur-workflow';

    FOR a_info IN
        SELECT * FROM qiita.artifact WHERE command_id = old_deblur_id
    LOOP
        SELECT processing_job_id INTO pj_id
            FROM qiita.processing_job
                JOIN qiita.artifact_output_processing_job USING (processing_job_id)
            WHERE artifact_id = a_info.artifact_id;

        UPDATE qiita.processing_job
            SET command_id = new_deblur_id
            WHERE processing_job_id = pj_id;

        UPDATE qiita.artifact
            SET command_id = new_deblur_id
            WHERE artifact_id = a_info.artifact_id;

    END LOOP;
END $do$;

-- after initial deploy we realized that we need to make sure that all
-- pos-ref-db-fp  /databases/gg/13_8/sortmerna/88_otus
-- pos-ref-fp     /databases/gg/13_8/rep_set/88_otus.fasta

DO $do$
DECLARE
    deblur_id     BIGINT;
    a_info        RECORD;
    pj_id         UUID;
    parameters    JSON;
BEGIN
    SELECT command_id INTO deblur_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'deblur' AND s.version = '1.0.3' AND sc.name = 'deblur-workflow';

    FOR a_info IN
        SELECT * FROM qiita.artifact WHERE command_id = deblur_id
    LOOP
        parameters := ('{"Jobs to start": 5, "Indexed positive filtering database": "/databases/gg/13_8/sortmerna/88_otus", ' ||
                        '"Threads per sample": 1, "Insertion/deletion (indel) probability": 0.01, "Negative filtering database": "default", ' ||
                        '"Mean per nucleotide error rate": 0.005, ' ||
                        '"Error probabilities for each Hamming distance": "1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005", ' ||
                        '"Indexed negative filtering database": "default", "Sequence trim length (-1 for no trimming)": -1, ' ||
                        '"Positive filtering database": "/databases/gg/13_8/rep_set/88_otus.fasta", "Maximum number of insertion/deletion (indel)": 3, ' ||
                        '"Minimum dataset-wide read threshold": 0, "Minimum per-sample read threshold": 2, "Demultiplexed sequences": "' ||
                        (a_info.command_parameters->>'seqs-fp')::varchar || '"}')::json;

        SELECT processing_job_id INTO pj_id
            FROM qiita.processing_job
                JOIN qiita.artifact_output_processing_job USING (processing_job_id)
            WHERE artifact_id = a_info.artifact_id;

        UPDATE qiita.processing_job
            SET command_parameters = parameters
            WHERE processing_job_id = pj_id;

        UPDATE qiita.artifact
            SET command_parameters = parameters
            WHERE artifact_id = a_info.artifact_id;

    END LOOP;
    UPDATE qiita.default_parameter_set
        SET parameter_set = ('{"Jobs to start": 5, "Indexed positive filtering database": "/databases/gg/13_8/sortmerna/88_otus", ' ||
                        '"Threads per sample": 1, "Insertion/deletion (indel) probability": 0.01, "Negative filtering database": "default", ' ||
                        '"Mean per nucleotide error rate": 0.005, ' ||
                        '"Error probabilities for each Hamming distance": "1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005", ' ||
                        '"Indexed negative filtering database": "default", "Sequence trim length (-1 for no trimming)": -1, ' ||
                        '"Positive filtering database": "/databases/gg/13_8/rep_set/88_otus.fasta", "Maximum number of insertion/deletion (indel)": 3, ' ||
                        '"Minimum dataset-wide read threshold": 0, "Minimum per-sample read threshold": 2}')::json
        WHERE command_id = deblur_id;
END $do$;
