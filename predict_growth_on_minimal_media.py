from __future__ import print_function
import sys
import pickle
from os import listdir
sys.path.insert(0, "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/")
import PyFBA


# Load the model seed database and change the incorrect reactions
compounds, reactions, enzymes =\
           PyFBA.parse.model_seed.compounds_reactions_enzymes('gramnegative')


# Load the set of reactions from original draft model
original_rxns = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive'
                             '/citrobacter_gapfilling_4/citrobacter_draft_reactions.p',
                             'rb'))

# Load the set of reactions added from gap-filling on LB media
gf_LB = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                         'citrobacter_gapfilling_4/ArgonneLB_added_reactions.p','rb'))

# Load the set of reactions added from gapfilling on minimal media
"""
gf_min_media = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                                'citrobacter_gapfilling_4/min_media_gapfill_added_rxns.p',
                                'rb'))
"""
gf50_percent = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_more_than_half_media.p',
                    'rb'))

gf25_percent = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_more_than_quarter_media.p',
                    'rb'))

gf10_percent = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_more_than_10_percent_media.p',
                    'rb'))

gf5_percent = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_more_than_5_percent_media.p',
                    'rb'))

gf_2_media = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_only_2_media.p',
                    'rb'))

gf_1_media = pickle.load(open('/Users/Taylor/anthill_backup/backup_archive/'
                    'citrobacter_gapfilling_4/reactions_added_on_only_1_media.p',
                    'rb'))

# For reactions added on only 1 or 2 media, only included the transport reactions
gf_rest_transport = set()
gf_rest_non_transport = set()
"""
for r in gf_2_media:
    if reactions[r].is_transport:
        gf_rest_transport.add(r)
    else:
        gf_rest_non_transport.add(r)
"""
for r in gf_1_media:
    if reactions[r].is_transport:
        gf_rest_transport.add(r)
    else:
        gf_rest_non_transport.add(r)

# Set the reactions to run in FBA
reactions_to_run = set()
reactions_to_run.update(original_rxns)
reactions_to_run.update(gf_LB)
#reactions_to_run.update(gf_min_media)
reactions_to_run.update(gf50_percent)
reactions_to_run.update(gf25_percent)
reactions_to_run.update(gf10_percent)
reactions_to_run.update(gf5_percent)
reactions_to_run.update(gf_2_media)
#reactions_to_run.update(gf_1_media)
reactions_to_run.update(gf_rest_transport)
#reactions_to_run.update(gf_rest_non_transport)

# Load the media conditions and experimental phenotypic growth data
exp_growth_results = {}
with open('/Users/Taylor/Desktop/c.sedlakii_growth.txt', 'r') as fin:
    for i, line in enumerate(fin):
        if i==0:
            continue
        condition, result = line.strip().split('\t')
        exp_growth_results[condition] = int(result)


# Set the biomass equation for FBA
biomass_equation = PyFBA.metabolism.biomass_equation('gramnegative')

        
# Dictionary to record FBA results on the various media
fba_growth_results = {}
# Iterate through all media conditions and run FBA
for media_condition in exp_growth_results:
    print("\nMEDIA CONDITION: " + media_condition)
    # Load the media for FBA
    media = PyFBA.parse.read_media_file(
        '/Users/Taylor/gapfilling_metabolic_networks/PyFBA/'
        'media/' + media_condition + '.txt')
    # Run the FBA
    status, value, growth =\
            PyFBA.fba.run_fba(compounds, reactions, reactions_to_run, media,
                          biomass_equation, verbose=True)
    print("FBA run on {} media has a biomass flux value"
      " of {} --> Growth: {}".format(media_condition, value, growth))
    
    # Record the FBA results
    fba_growth_results[media_condition] = int(growth)
    
"""
# Write FBA growth results to file
with open('citrobacter_gapfilling_2/c_sedlakii_initial_fba_growth_results.txt', 'w') as fout:
    fout.write('MEDIA\tGROWTH\n')
    for media in fba_growth_results:
        fout.write(media + '\t' + str(fba_growth_results[media]) + '\n')
"""

# Check if FBA results agree with experimental phenotypic growth data
results = {'tp': [], 'tn': [], 'fp': [], 'fn': []}
count_agree = 0
for media in fba_growth_results:
    # Results in agrement (tp and tn)
    if fba_growth_results[media] == exp_growth_results[media]:
        count_agree += 1
        if exp_growth_results[media] == 1:
            results['tp'].append(media)
        else:
            results['tn'].append(media)
    # Results that don't agree (fp and fn)
    else:
        if exp_growth_results[media] == 0:
            results['fp'].append(media)
        else:
            results['fn'].append(media)


percent_agreement = (float(count_agree) / len(fba_growth_results)) * 100
print('\n\nThe agreement between the FBA growth results and the experimental '
      'phenotypic growth results is {} %.'.format(percent_agreement))
print('\nTotal number of predictions: {}'.format(len(fba_growth_results)))
print('\nCOREECT PREDICTIONS:')
print('\tTP: {}'.format(len(results['tp'])))
print('\tTN: {}'.format(len(results['tn'])))
print('\nINCORRECT PREDICTIONS:')
print('\tFP: {}'.format(len(results['fp'])))
print('\tFN: {}'.format(len(results['fn'])))

print('Gap-filling is needed on {} media.'.format(len(results['fn'])))

# Save the dictionary of simulated growth results
#pickle.dump(results, open('citrobacter_gapfilling_2/c_sedlakii_initial_fba_growth_results.p','wb'))

    
    
    


    

    

    

    

    
