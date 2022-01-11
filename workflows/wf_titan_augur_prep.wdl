version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/task_versioning.wdl" as versioning

# goal: add more metadata fields.

workflow titan_augur_prep {
    input {
        String    assembly
        String    collection_date
        String    iso_country
        String    iso_state
        String    iso_continent
        String?   iso_county
        String    pango_lineage
        String?   extra_field1
        String?   extra_field1_title
        String?   extra_field2
        String?   extra_field2_title
        String?   extra_field3
        String?   extra_field3_title
        

    }

    call nextstrain.prep_augur_metadata {
        input:
            assembly=assembly,
            collection_date = collection_date,
            iso_country = iso_country,
            iso_state = iso_state,
            iso_continent = iso_continent,
            iso_county = iso_county,
            pango_lineage = pango_lineage,
            extra_field1 = extra_field1,
            extra_field2 = extra_field2,
            extra_field3 = extra_field3,
            extra_field1_title = extra_field1_title,
            extra_field2_title = extra_field2_title,
            extra_field3_title = extra_field3_title
            
        
    }
    call versioning.version_capture{
      input:
    }
    output {
      String titan_augur_run_version            = version_capture.phvg_version
      String titan_augur_run_analysis_date      = version_capture.date
      File   augur_metadata                     = prep_augur_metadata.augur_metadata

    }
}
