from __future__ import print_function
import pickle
from likelihood_gapfill import suggest_additional_reactions
import sys
sys.path.insert(0, "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/")
import PyFBA


# Load the Model SEED database
compounds, reactions, enzymes =\
           PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')


# Load the draft model reactions and roles
draft_rxns = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                              'citrobacter_gapfilling_4/citrobacter_draft_reactions.p',
                              'rb'))
draft_roles = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                              'citrobacter_gapfilling_4/citrobacter_draft_roles.p',
                              'rb'))
# Load the reactions and roles added in gap-filling on rich media and add them to
# the draft reactions and draft roles
LB_added_rxns = pickle.load(
    open('/Users/Taylor/anthill_backup/backup_archive/citrobacter_gapfilling_4/'
        'ArgonneLB_added_reactions.p','rb'))
draft_rxns.update(LB_added_rxns)

LB_added_roles = set()
rxn2role = PyFBA.filters.reactions_to_roles(LB_added_rxns)
for i in rxn2role:
    LB_added_roles.update(rxn2role[i])
draft_roles.update(LB_added_roles)


# Read in media conditions in which the organism is known to grow
pos_growth_media = set()
with open('/Users/Taylor/Desktop/citrobacter_sedlakii/'
          'c.sedlakii_pos_growth.txt','r') as fin:
    for line in fin:
        pos_growth_media.add(line.strip())
print("{} growth media conditions to gap-fill on".format(len(pos_growth_media)))


#Set the biomass equation
biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative')

reactions_suggested_per_media = {}
all_suggested_rxns = set()
all_suggested_roles = set()
all_suggested_rxn_source = {}
# Suggest reactions and functional roles possibly missing form the model for
# each of the minimal media sources that gap-filling will be performed on
for media_condition in pos_growth_media:
    print("\n\n\nSuggesting reactions on {} media...".format(media_condition))

    # Read the media file and set the media variable
    media = PyFBA.parse.read_media_file\
            ('/Users/Taylor/gapfilling_metabolic_networks/'
             'PyFBA/media/' + media_condition + '.txt')

    # Suggest additional reactions
    suggested_rxns, suggested_roles, source = suggest_additional_reactions(compounds,
            reactions, draft_rxns, draft_roles, media, biomass_equation,
            "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/example_data/"
            "Citrobacter/ungapfilled_model/closest.genomes.roles",
            "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/example_data/"
            "Citrobacter/ungapfilled_model/citrobacter.roles")

    # Record which reactions were added from the media
    reactions_suggested_per_media[media_condition] = suggested_rxns

    # Update the set of suggested reactions
    all_suggested_rxns.update(suggested_rxns)
    # Record the source of the suggested reaction
    for rxn in suggested_rxns:
        if rxn not in all_suggested_rxn_source:
            all_suggested_rxn_source[rxn] = [source[rxn]]
        else:
            all_suggested_rxn_source[rxn].append(source[rxn])

    # Update the set of suggested roles
    all_suggested_roles.update(suggested_roles)



print("{} reactions were suggested in total from all of the minimal"
      " media sources.".format(len(all_suggested_rxns)))
print("{} roles were suggested in total from all of the minimal media sources."
      .format(len(all_suggested_roles)))

# Save the suggested reactions and roles
pickle.dump(all_suggested_rxns, open('citrobacter_gapfilling_4/all_suggested_reactions.p','wb'))
pickle.dump(all_suggested_roles, open('citrobacter_gapfilling_4/all_suggested_roles.p','wb'))
pickle.dump(all_suggested_rxn_source, open('citrobacter_gapfilling_4/all_suggested_reactions_source.p','wb'))
pickle.dump(reactions_suggested_per_media, open('reactions_suggested_per_min_media.p','wb'))









