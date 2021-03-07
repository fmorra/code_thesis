def region_counter(region_list, substring):

    import numpy as np

    # substrings_found = np.zeros((len(region_list),), dtype="int")
    substrings_found = [i for i in region_list if substring in i]
    # for string_ind in range(len(region_list)):
    #     find_result = region_list[string_ind].find(substring)
    #     substrings_found[string_ind] = find_result
    #     substrings_found[substrings_found > -1] = 1
    #     substrings_found[substrings_found == -1] = 0

    return substrings_found
