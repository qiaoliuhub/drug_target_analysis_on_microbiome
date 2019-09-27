#! /usr/bin/env python
import config
import pandas as pd
import logging

# Setting up log file
formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s %(name)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S')
fh = logging.FileHandler(config.logfile, mode='a')
fh.setFormatter(fmt=formatter)
logger = logging.getLogger("Select Drugs")
logger.addHandler(fh)
logger.setLevel(logging.DEBUG)

### drugs inhibit "bad"
### drugs inhibit "good"
### drugs activate "bad"
### drugs activate "good"
if __name__ == "__main__":

    beneficial_species = pd.read_csv(config.beneficial_species_file)
    pathogen_species = pd.read_csv(config.pathogen_species_file)
    logger.debug("Successfully read in species")
    activation = pd.read_csv(config.activation, header = None)
    logger.debug("Successfully read in activation")
    inhibition = pd.read_csv(config.inhibition, header = None)
    logger.debug("Successfully read in inhibition")
    specie_name_id = pd.read_csv(config.specie_id_map)

    pathogen_species_id = set(pathogen_species['id'])
    pathogen_taxon_id = set(specie_name_id[specie_name_id['id'].isin(pathogen_species_id)]['taxon_id'])
    beneficial_species_id = set(beneficial_species['id'])
    beneficial_taxon_id = set(specie_name_id[specie_name_id['id'].isin(beneficial_species_id)]['taxon_id'])

    path_acti = activation[activation[6].isin(pathogen_taxon_id)]
    bene_acti = activation[activation[6].isin(beneficial_taxon_id)]
    path_inhi = inhibition[inhibition[6].isin(pathogen_taxon_id)]
    bene_inhi = inhibition[inhibition[6].isin(beneficial_taxon_id)]

    output_file = iter([config.path_acti_bene_acti, config.path_acti_bene_inhi, config.path_inhi_bene_acti, config.path_inhi_bene_inhi])
    path_affe, bene_affe = set(), set()
    for path in [path_acti, path_inhi]:
        for bene in [bene_acti, bene_inhi]:
            path_chem = set(x for x in (set(path[0]) | set(path[1])) if x.startswith('CID'))
            path_affe |= path_chem
            bene_chem = set(x for x in (set(bene[0]) | set(bene[1])) if x.startswith('CID'))
            bene_affe |= bene_chem
            pd.DataFrame({'CIDs': list(path_chem - bene_chem)}).to_csv(next(output_file), index = False)
            logger.debug("Successfully write out files")

    pd.DataFrame({'CIDs': list(path_affe-bene_affe)}).to_csv(config.path_affe_bene_affe, index = False)
