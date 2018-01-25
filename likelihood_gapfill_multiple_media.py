from __future__ import print_function
import pickle
from likelihood_gapfill import build_draft_model, suggest_additional_reactions, likelihood_gapfill_optimization
import sys
sys.path.insert(0, "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/")
import PyFBA

# Load the Model SEED database
compounds, reactions, enzymes =\
           PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')

# Get the set of essential reactions that are present in all models
essentials = PyFBA.gapfill.suggest_essential_reactions()

"""
# Build the draft model for the organism
draft_roles, draft_rxns =\
            build_draft_model('/Users/Taylor/Desktop/citrobacter_sedlakii/'
                              '67826.8/assigned_functions.txt')
"""

# Load the draft model reactions and roles
draft_rxns = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                              'citrobacter_gapfilling_4/citrobacter_draft_reactions.p',
                              'rb'))
draft_roles = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                              'citrobacter_gapfilling_4/citrobacter_draft_roles.p',
                              'rb'))
# Load the reactions and roles added in gap-filling on rich media and add them to
# the draft reactions and draft roles
LB_added_rxns =\
    pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/ArgonneLB_added_reactions.p','rb'))
draft_rxns.update(LB_added_rxns)

LB_added_roles = set()
rxn2role = PyFBA.filters.reactions_to_roles(LB_added_rxns)
for i in rxn2role:
    LB_added_roles.update(rxn2role[i])
draft_roles.update(LB_added_roles)


#Set the biomass equation
biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative')

# Read in media conditions in which the organism is known to grow
pos_growth_media = set()
with open('/Users/Taylor/Desktop/citrobacter_sedlakii/'
          'c.sedlakii_pos_growth.txt','r') as fin:
    for line in fin:
        pos_growth_media.add(line.strip())
print("{} growth media conditions to gap-fill on".format(len(pos_growth_media)))


gapfill_added_rxns = set()
gapfill_media_source = {}
# Run gap-filling on each of the media conditions
for media_condition in pos_growth_media:
    print("\n\n\nGap-filling on {} media...".format(media_condition))

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
    print("\n{} reactions were suggested to complete the model for {} media.\n"
          .format(len(suggested_rxns), media_condition))
                               
    # Run likelihood gap-filling optimization
    added_reactions, added_rxn_fluxes =\
        likelihood_gapfill_optimization(compounds, reactions,
             draft_rxns, suggested_rxns, biomass_equation,
             media, "/Users/Taylor/anthill_backup/backup_archive/"
             "genome_reaction_probabilities.txt", essentials)
    print("\n{} reactions were added in gap-filling on {} media.\n"
          .format(len(added_reactions), media_condition))

    # Write out the gapfill reactions added to the model on the media
    # condition to a text file
    fout = open("citrobacter_gapfilling_4/min_media_gapfilling_solutions/"
                "gapfill_reactions_" + media_condition + ".txt", "w")
    for rxn in added_reactions:
        fout.write(rxn + "\n")
    fout.close()
    
    # Record reactions added to the model in gapfilling
    gapfill_added_rxns.update(added_reactions)
    # Record which media the reaction was added from
    for rxn in added_reactions:
        if rxn not in gapfill_media_source:
            gapfill_media_source[rxn] = [media_condition]
        else:
            gapfill_media_source[rxn].append(media_condition)


# Save all of the reactions added in gapfilling on the multiple media types
print("\n\n\n{} reactions in total were added in gap-filling on the {} media "
      "conditions.".format(len(gapfill_added_rxns), len(pos_growth_media)))
pickle.dump(gapfill_added_rxns, open("citrobacter_gapfilling_4/"
                                     "min_media_gapfill_added_rxns.p","wb"))
pickle.dump(gapfill_media_source, open("citrobacter_gapfilling_4/"
                                       "min_media_added_reaction_media_source.p",
                                       "wb"))

