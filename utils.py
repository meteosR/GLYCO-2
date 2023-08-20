
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
        
    elif submission_data["module_type"] == "subset":
        if submission_data["residue_list_file"] is None:
            valid_submission = False
            error_message_validation = "Subset module requires a residue list file."
    else:
        valid_submission = False
        error_message_validation = "Incorrect module type."

    if submission_data["glycan_names"] is None or submission_data["glycan_names"].strip() == "":
        valid_submission = False
        error_message_validation = "Incorrect glycan names."
    else:
        submission_data["glycan_names"] = [x.strip() for x in submission_data["glycan_names"].split(",")]

    if float(submission_data["cylinder_radius"]) <= 0:
        valid_submission = False
        error_message_validation = "Cylinder radius must be > 0"

    return valid_submission, error_message_validation
