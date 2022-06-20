import random
def structure_catalogue_maker(Number_of_structures, 
                              Number_of_atoms, 
                              lower_atom_number, 
                              higher_atom_number,
                             ) -> list:
    """
    This function makes a shuffled list containing 'Number_of_structures' number of lists which each is 
    'Number_of_atoms' long and is randomly distributed with 0's and 1's whereas the minimum number of 1's are 
    'lower_atom_number' and the maximum number of 1's are 'higher_atom_number'.
    
    returns catalogue as list
    """
    
    print (f"Starting to make a structure catalogue with: {Number_of_structures} structure from the starting model.")
    print (f"The structure will have between {lower_atom_number} and {higher_atom_number} atoms")
    
    # generate empty list
    structure_catalogue = []
    
    for i in range(Number_of_structures):
        one_count = random.randint(lower_atom_number, higher_atom_number)
        zero_count = Number_of_atoms  - one_count
        
        # create a list of 0 and 1's
        my_list = [0]*zero_count + [1]*one_count
        
        # shuffle the list
        random.shuffle(my_list)
        my_list.insert(0, one_count)
        
        # append to catalogue
        structure_catalogue.append(my_list)
        
    print("Permutations Succeeded")
    return structure_catalogue