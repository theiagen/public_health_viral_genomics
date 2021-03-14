version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow titan_augur_prep {
	input {
		String 		assembly
		String    collection_date
		String    iso_country
		String    iso_state
		String    iso_continent
		String?   iso_county
		String    pangolin_lineage

	}

	call nextstrain.prep_augur_metadata {
		input:
			assembly=assembly,
			collection_date = collection_date,
			iso_country = iso_country,
			iso_state = iso_state,
			iso_continent = iso_continent,
			iso_county = iso_county,
			pangolin_lineage = pangolin_lineage
	}

	output {
			File     augur_metadata = prep_augur_metadata.augur_metadata

	}
}
