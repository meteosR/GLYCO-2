
def data_validator(submission_data):
    error_message_validation = ""
    valid_submission = True
    
    # Validate protein parameters #
    
    if submission_data["module_type"] == "all_atom":
        # Make sure BSA threshold present
        if type(submission_data["surface_threshold"]) is None:
            valid_submission = False
            error_message_validation = "Minimum RSA is not valid"
        else:  
            if float(submission_data["surface_threshold"]) < 0:
                valid_submission = False
                error_message_validation = "Minimum RSA must be >= 0"
                
            if float(submission_data["probe_radius"]) <= 0:
                valid_submission = False
                error_message_validation = "Probe radius must be > 0"
        
        # Check selection method (chain or type)
        if submission_data["selection_type_p"] == "chain":
            
            if submission_data["protein_chains"] is None or submission_data["protein_chains"].strip() == "":
                valid_submission = False
                error_message_validation = "Incorrect protein chains."
            else: 
                submission_data["protein_chains"] = submission_data["protein_chains"].split(",")
            
        elif submission_data["selection_type_p"] == "type":  
            # No further input required, will use default protein alphabet      
            pass         
        else:
            valid_submission = False
            error_message_validation = "Incorrect selection type for protein."
        
    elif submission_data["module_type"] == "subset":
        if submission_data["residue_list_file"] is None:
            valid_submission = False
            error_message_validation = "Subset module requires a residue list file."
    else:
        valid_submission = False
        error_message_validation = "Incorrect module type."
                
    # Validate glycan paramters #
    if submission_data["selection_type_g"] == "chain":
        if submission_data["glycan_chains"] is None or submission_data["glycan_chains"].strip() == "":
            valid_submission = False
            error_message_validation = "Incorrect glycan chains."
        else:
            submission_data["glycan_chains"] = [x.strip() for x in submission_data["glycan_chains"].split(",")]
            
    elif submission_data["selection_type_g"] == "type":  
        if submission_data["glycan_names"] is None or submission_data["glycan_names"].strip() == "":
            valid_submission = False
            error_message_validation = "Incorrect glycan names."
        else:
            submission_data["glycan_names"] = [x.strip() for x in submission_data["glycan_names"].split(",")]   
                 
    else:
        valid_submission = False
        error_message_validation = "Incorrect selection type for glycans."    
        
        
    if float(submission_data["cylinder_radius"]) <= 0:
        valid_submission = False
        error_message_validation = "Cylinder radius must be > 0"
        
    

    return valid_submission, error_message_validation