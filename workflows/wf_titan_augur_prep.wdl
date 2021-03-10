version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow mercury_se_prep {
	input {
		String    submission_id
		String    collection_date
		String    iso_country
		String    iso_state
		String    iso_continent
		String    pangolin_lineage

	}

	call nextstrain.prep_augur_metadata {
		input:
			submission_id = submission_id,
			collection_date = collection_date,
			iso_country = iso_country,
			iso_state = iso_state,
			iso_continent = iso_continent,
			pangolin_lineage = pangolin_lineage
	}

	output {
			File     augur_metadata = "${submission_id}.augur_metadata.csv"

	}
}
