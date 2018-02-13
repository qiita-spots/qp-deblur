-- moving deblur v1.0.3 artifacts to v1.0.4
-- we need to add: 'Reference phylogeny for SEPP': 'Greengenes_13.8'

DO $do$
DECLARE
    old_deblur_id BIGINT;
    new_deblur_id BIGINT;
    a_info        RECORD;
    pj_id         UUID;
    parameters    TEXT;
    params        JSON;
BEGIN
    SELECT command_id INTO old_deblur_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'deblur' AND s.version = '1.0.3' AND sc.name = 'deblur';

    SELECT command_id INTO new_deblur_id
        FROM qiita.software_command sc
            JOIN qiita.software s USING (software_id)
        WHERE s.name = 'deblur' AND s.version = '1.0.4' AND sc.name = 'deblur';


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

        SELECT processing_job_id INTO pj_id
            FROM qiita.processing_job
                JOIN qiita.artifact_output_processing_job USING (processing_job_id)
            WHERE artifact_id = a_info.artifact_id;

        parameters := LEFT((a_info.command_parameters)::VARCHAR, -1);
        params := (parameters || ', "Reference phylogeny for SEPP": "Greengenes_13.8"}')::JSON;

        UPDATE qiita.processing_job
            SET command_parameters = params
            WHERE processing_job_id = pj_id;

        UPDATE qiita.artifact
            SET command_parameters = params
            WHERE artifact_id = a_info.artifact_id;

    END LOOP;
    UPDATE qiita.command_parameter
      SET check_biom_merge = True
      WHERE parameter_name = 'Reference phylogeny for SEPP'
        AND command_id=new_deblur_id;
END $do$;
