from __future__ import print_function
import copy
import sys
sys.path.insert(0, "/Users/Taylor/gapfilling_metabolic_networks/PyFBA/")
import PyFBA




def build_draft_model(assigned_functions_file, orgtype='gramnegative', verbose=True):
    """
    Bild a draft metabolic model from the RAST assigned functions
    file.

    Read in the functional roles from the assigned functions file,
    map the roles to the associated reactions using hte mappings
    from the Model SEED database, and return the set of
    functional roles and the set of reactions for the organism.

    :param assigned_functions_file: Filepath to the RAST assigned functions file
    :type assigned_functions_file: string
    :param orgtype: Organism type 
    :type orgtype: string
    :param verbose: Verbose output
    :type verbose: bool
    :return: A set of functional roles and a set of reaction ids
    :rtype: (set, set)
    """
    
    # Read RAST assigned functions file to get functional roles in the genome
    assigned_fxns = PyFBA.parse.read_assigned_functions(assigned_functions_file)
    if verbose:
        print("Loading initial draft model for organism ...")
    roles_present = set()
    for roles in assigned_fxns.values():
        roles_present.update(roles)
    if verbose:
        print("There are {} unique roles in this genome".format(len(roles_present)))

    # Convert functional roles into reactions
    roles_to_reactions = PyFBA.filters.roles_to_reactions(roles_present)
    reactions_present = set()
    for role in roles_to_reactions:
        reactions_present.update(roles_to_reactions[role])
    if verbose:
        print("There are {}".format(len(reactions_present)),
              "unique reactions associated with this genome.")

    # Read in ModelSEED database of compounds, reactions, and enzyme complexes
    compounds, reactions, enzymes =\
            PyFBA.parse.model_seed.compounds_reactions_enzymes(orgtype)

    # Update the reactions to run, making sure that all the reactions
    # are in our reactions database
    tempset = set()
    for r in reactions_present:
        if r in reactions:
            tempset.add(r)
        else:
            if verbose:
                print("Reaction ID {}".format(r),
                      "is not in our reactions list. Skipped", file=sys.stderr)
    reactions_present = tempset

    return roles_present, reactions_present




def suggest_additional_reactions(compounds, reactions, draft_reactions,
                                  draft_roles, media, biomass_equation,
                                  close_roles_file, genus_roles_file,
                                  verbose=True):
    """
    Suggest additional reactions to add to a draft model to enable the model
    to grow on a media type where it is known to grow.  Reactions are suggested
    using functions from the PyFBA gapfill module.  Reactions are suggested
    from the following sources:
        - Reactions based on the media
        - Reactions from roles present in RAST close genomes
        - Reactions from roles present in genomes from the same genus
        - Essential reactions that are present in all models
        - Reactions that complete subsystems
        - Reactions connected to orphan compounds
        - Reactions that connect to other compounds present in the model
          and meet a chosen threshold for probability of existence (probability
          is based upon the fraction of compounds for the reaction that are
          present in the model)

    :param compounds: The dictionary of all compounds from Model SEED
    :type compounds: dict
    :param reactions: The dictionary of all reactions from Model SEED
    :type reactions: dict
    :param draft_reactions: A set of reaction ids for the reactions present in the draft model
    :type draft_reactions: set
    :param draft_roles: A set of functional roles present in the draft model
    :type draft_roles: set
    :param media: A set of compounds present in the media
    :type media: set
    :param biomass_equation: The biomass equation as a Reaction object
    :type biomass_equation: metabolism.Reaction object
    :param close_roles_file: A filepath to a file with a list of roles present in RAST close genomes
    :type close_roles_file: string
    :param genus_roles_file: A filepath to a file with a list of roles present in genomes from the same genus
    :type genus_roles_file: string
    :param verbose: Verbose output
    :type verbose: bool
    :return: A set of reactions possibly missing from the model, a set of roles possibly missing
        from the model, and a dictionary of source for the missing reactions
    :rtype: (set, set, dict)
    """

    # Initialize the reactions to run as the set of reactions from the draft model
    reactions_to_run = copy.copy(draft_reactions)
    
    # TEST IF DRAFT MODEL GROWS ON THE MEDIA
    status, value, growth =\
            PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                 media, biomass_equation)
    print("Initial FBA run has a biomass flux value"
          " of {} --> Growth: {}".format(value, growth))


    # PROPOSE REACTIONS TO ADD TO THE MODEL IF IT WON'T GROW ON THE MEDIA
    if not growth:
        # Keep track of the suggested reactions
        added_reactions = []
        reaction_source = {}
        
        # SUGGEST REACTIONS FROM MEDIA
        if verbose:
            print("\nFinding reactions from media...")
        media_reactions = PyFBA.gapfill.suggest_from_media(compounds, reactions,
                                                            reactions_to_run, media)
        added_reactions.append(("media", media_reactions))
        reactions_to_run.update(media_reactions)
        for rxn in media_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'media_reactions'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding media reactions, the biomass reaction "
                  "has a flux of {} --> Growth: {}".format(value, growth))
    

    if not growth:
        # SUGGEST REACTIONS FROM RAST CLOSELY-RELATED ORGANISMS
        if verbose:
            print("\nFinding reactions from closely-related organisms...")
        close_reactions = PyFBA.gapfill.suggest_from_roles(close_roles_file, reactions)
        # Find which of the suggested reactions are new
        close_reactions.difference_update(reactions_to_run)
        added_reactions.append(("close genomes", close_reactions))
        reactions_to_run.update(close_reactions)
        for rxn in close_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'close_genomes'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding reactions from RAST close genomes, "
                  "the biomass reaction has a flux of {} --> Growth: {}".format(value, growth))

    if not growth:
        # SUGGEST REACTIONS FROM ALL GENOMES IN SAME GENUS
        if verbose:
            print("\nFinding reactions from species in same genera...")
        genus_reactions = set()
        genus_reactions = PyFBA.gapfill.suggest_from_roles(genus_roles_file, reactions)
        # Find which of the suggested reactions are new
        genus_reactions.difference_update(reactions_to_run)
        added_reactions.append(("genus_reactions", genus_reactions))
        reactions_to_run.update(genus_reactions)
        for rxn in genus_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'genus_reactions'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding reactions from other species in the same genus, "
                  "the biomass reaction has a flux of {} --> Growth: {}".format(value, growth))

    if not growth:
        # SUGGEST ESSENTIAL REACTIONS
        if verbose:
            print("\nFinding essential reactions...")
        essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
        # Find which of the suggested reactions are new
        essential_reactions.difference_update(reactions_to_run)
        added_reactions.append(("essential", essential_reactions))
        reactions_to_run.update(essential_reactions)
        for rxn in essential_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'essential_ractions'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding essential reactions, the biomass reaction has"
                  " a flux of {} --> Growth: {}".format(value, growth))

    if not growth:
        # SUGGEST REACTIONS THAT COMPLETE SUBSYSTEMS
        if verbose:
            print("\nFinding reactions that complete subsystems...")
        subsystem_reactions =\
                            PyFBA.gapfill.suggest_reactions_from_subsystems(reactions,
                                                                             reactions_to_run,
                                                                             threshold=0.5)
        added_reactions.append(("subsystems", subsystem_reactions))
        reactions_to_run.update(subsystem_reactions)
        for rxn in subsystem_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'subsystem_reactions'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding subsystem reactions, the biomass reaction "
                  "has a flux of {} --> Growth: {}".format(value, growth))
    
    if not growth:
        # SUGGEST REACTIONS THAT CONNECT TO ORPHAN COMPOUNDS
        if verbose:
            print("\nFinding reactions conecting to orphan compounds...")
        orphan_reactions = PyFBA.gapfill.suggest_by_compound(compounds, reactions,
                                                              reactions_to_run,
                                                              max_reactions=1)
        added_reactions.append(("orphans", orphan_reactions))
        reactions_to_run.update(orphan_reactions)
        for rxn in orphan_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'orphan_compounds'

        # Test for growth
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding reactions connecting to orphan compounds, "
                  "the biomass reaction has a flux of {} --> Growth: {}".format(value, growth))

    if not growth:
        # SUGGEST COMPOUND-PROBABILITY REACTIONS
        if verbose:
            print("\nFinding compound-probability reactions...")
        probable_reactions = PyFBA.gapfill.compound_probability(reactions,
                                                                 reactions_to_run,
                                                                 cutoff=0,
                                                                 rxn_with_proteins=True)
        probable_reactions.difference_update(reactions_to_run)
        added_reactions.append(("compound probability", probable_reactions))
        reactions_to_run.update(probable_reactions)
        for rxn in probable_reactions:
            if rxn not in reaction_source:
                reaction_source[rxn] = 'probable_reactions'
    
        status, value, growth =\
                PyFBA.fba.run_fba(compounds, reactions, reactions_to_run,
                                     media, biomass_equation)
        if verbose:
            print("After adding reactions based on compound probability, "
                  "the biomass reaction has a flux of {} --> Growth: {}".format(value, growth))


    # GET THE SET OF REACTIONS THAT MAY NEED TO BE ADDED TO THE MODEL
    missing_reactions = set()
    for i in added_reactions:
        missing_reactions.update(i[1])
    if verbose:
        print("\nThere are {} reactions that may need to be added to the model"
              .format(len(missing_reactions)))

    # MAP THE REACTIONS THAT MAY NEED TO BE ADDED TO THEIR ASSOCIATED
    # FUNCTIONAL ROLES
    missing_reactions_to_roles = PyFBA.filters.reactions_to_roles(missing_reactions)
    missing_roles = set()
    for r in missing_reactions_to_roles:
        missing_roles.update(missing_reactions_to_roles[r])
    if verbose:
        print("\nThere are {} functional roles that may be missing from the model."
              .format(len(missing_roles)))

    return missing_reactions, missing_roles, reaction_source




def likelihood_gapfill_optimization(compounds, reactions, original_reactions,
                                     suggested_reactions, biomass_equation,
                                     media, role_probabilities_file,
                                     essential_reactions, verbose=True):
    """
    Run FBA in the likelihood-based gapfill mode to determine which of the
    suggested reactions to add to the draft model to enable growth. Reactions
    with higher associated probabilities are favored for addition to the model
    due to the way objective coefficients are calculated fromt the reaction
    probabilites and the objective function utilized in the optimization.

    :param compounds: The dictionary of compounds from the Model SEED database
    :type compounds: dict
    :param reactions: The dictionary of reactions from the Model SEED database
    :type reactions: dict
    :param original_reactions: The set of reaction ids from the draft model
    :type original_reactions: set
    :param suggested_reactions: The set of suggested reaction ids to be used in
        gapfilling to complete the model and enable growth
    :type suggested_reactions: set
    :param biomass_equation: The biomass equation as a Reaction object
    :type biomass_equation: metabolism.Reaction object
    :param media: A set of compounds present in the media
    :type media: set
    :param role_probabilities_file: Filepath to file containing a list of roles
        and associated probabilities for the roles
    :type role_probabilities_file: string
    :param essential_reactions: The set of essntial reactions (returned by the PyFBA.gapfill.suggest_essential_reactions() function)
    :type essential_reactions: set
    :param verbose: Verbose output
    :type verbose: bool
    :return: A set of reactions added in the gapfilling process and a dictionary
        of fluxes for the added reactions
    :rtype: (set, dict)
    """

    # Read in the role probabilities from text file
    if verbose:
        print("Loading reaction probabilities file ...")
    rxn_probs = {}
    with open(role_probabilities_file,'r') as fin:
        for i, line in enumerate(fin):
            if i==0:
                continue
            r, p = line.strip().split('\t')
            rxn_probs[r] = float(p)

    
    # Enforce that all the candidate gap-filling reactions not present in the
    # original model can run only in the left to right (>) direction.  Reverse
    # all of the reacions that run right to left, and split the bidirectional
    # reactions into two separate left to right reactions.  The reactions
    # present in the original model can be left unmodified.
    if verbose:
        print("Enforcing all potential gapfilling reactions to run in the "
              " left to right direction ...")
    to_delete = []
    to_add = []
    rev_count = 0
    split_count = 0
    for rxn in suggested_reactions:
        # If the reaction is a transport reaction that runs left to right,
        # explicitly set the reaction bounds
        if reactions[rxn].direction == '>' and reactions[rxn].is_transport:
            reactions[rxn].lower_bound = 0.0
            reactions[rxn].upper_bound = 1000.0

        # Reverse the reaction if it runs right to left
        if reactions[rxn].direction == '<':
            rev_count += 1
            reactions[rxn].reverse_reaction()
            if reactions[rxn].is_transport:
                # Explicitly set the reaction bounds for transport reactions
                reactions[rxn].lower_bound = 0.0
                reactions[rxn].upper_bound = 1000.0

        # Split bidirectional reactions into two forward reactions:
        # one for left to right and one for right to left
        if reactions[rxn].direction == '=':
            split_count += 1
            # Split the reaction
            (fwd, rev) = reactions[rxn].split_reaction()
            if reactions[rxn].is_transport:
                # Explicitly set the reaction bounds for transport reactions
                fwd.lower_bound = 0.0
                fwd.upper_bound = 1000.0
                rev.lower_bound = 0.0
                rev.upper_bound = 1000.0
            # Add the forward and reverse reactions to the master
            # reactions dictionary
            reactions[fwd.name] = fwd
            reactions[rev.name] = rev
            # Keep track of reactions to delete from and add
            # to reactions_to_run set
            to_delete.append(rxn)
            to_add += [fwd.name, rev.name]

    if verbose:
        print("{} reactions reversed to run left to right.".format(rev_count))
        print("{} reactions split into separate forward and reverse reactions"
              .format(split_count))

   
    # Update the suggested_reactions set
    # Remove the bidirectional reactions from the suggested_reactions set
    for rxn in to_delete:
        suggested_reactions.discard(rxn)
    # Add the new forward and reverse reactions for the bidirectional
    # reactions that were split to the suggested_reactions set
    suggested_reactions.update(to_add)

    # Update the reactions probability hash with probabilities for the
    # forward and reverse reactions that were created from the
    # bidirectional reactions that were split
    for rxn in to_add:
        bidirect_rxn = rxn[:-2]
        if bidirect_rxn in rxn_probs:
            rxn_probs[rxn] = rxn_probs[bidirect_rxn]
    

    # Set the reactions to run in FBA
    reactions_to_run = set()
    reactions_to_run.update(original_reactions)
    reactions_to_run.update(suggested_reactions)

    # Run gapfilling by running the FBA in the likelihood gapfilling mode
    status, value, growth = \
            PyFBA.fba.run_fba(compounds, reactions, reactions_to_run, media,
                              biomass_equation, verbose=True,
                              likelihood_gapfill=True,
                              reaction_probs = rxn_probs,
                              original_reactions_to_run = original_reactions,
                              essential_reactions = essential_reactions)

    # Get the reaction fluxes from FBA to see which reactions ran
    reaction_flux = PyFBA.fba.reaction_fluxes()
    biomass_flux = reaction_flux["BIOMASS_EQN"]
    if biomass_flux >= 1.0:
        growth = True
    else:
        growth = False
    if verbose:
        print("After gap-filling, the biomass reaction has a flux "
              " of {} --> Growth: {}".format(biomass_flux, growth))

    # Check how many reactions ran and record which reactions were added
    # in gap-filling along with their fluxes
    rxns_running_count = 0
    gf_rxn_fluxes = {}
    for r in reaction_flux:
        if r == "BIOMASS_EQN":
            continue
        if reaction_flux[r] != 0.0 and "UPTAKE_SECRETION_REACTION" not in r:
            rxns_running_count += 1
            if r not in original_reactions:
                gf_rxn_fluxes[r.replace('_f','').replace('_r','')] = reaction_flux[r]
    # Record set of gapfilling reactions added
    gf_added_reactions = set(
        [i.replace('_f','').replace('_r','') for i in gf_rxn_fluxes.keys()])


    return gf_added_reactions, gf_rxn_fluxes

    

    




