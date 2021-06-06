def dat_processor(fname):
    # This function reads the .dat file containing node coordinates and the node labels making up each element and
    # stores the lines where they are defined for the next processing operations. Another .dat file, easier to read
    # for the next scripts, is also saved
    # Open the dat file
    with open(fname, "r+") as read_dat:
        # Read the file as a list of lines
        opened_dat = read_dat.readlines()
        print("Processing the dat file to make it more readable")
        # Initialize counters containing a series of lines in the dat file where we will perform the search for node
        # coordinates and element nodes
        node_lines = []
        elem_lines = []
        part_lines = []
        nset_lines = []
        node_keyword_counter = 0
        elem_keyword_counter = 0
        part_keyword_counter = 0
        nset_keyword_counter = 0
        file_identifiers = []
        # Define the keywords used to search for the lines delimiting the start and end of those regions
        node_search_key = '*Node'
        element_search_key = '*Element'
        part_search_key = 'PART INSTANCE'
        region_finisher_key = '*Nset'
        input_paragraph_end = 'OPTIONS BEING PROCESSED'
        line_counter = 0
        end_line_counter = 0
        # Iterate over each line of the dat file to find whether the keywords or keyphrases are contained in it, and if
        # so, save the number of those lines. Again, extract the names of all of the model parts as they will be used to
        # name the csv files containing the node coordinates and element nodes. The dat file only contains the part
        # values, without being divided into the different regions.
        for line in opened_dat:
            line = line.strip()
            line_counter += 1
            node_keyword_flag = line.find(node_search_key)
            elem_keyword_flag = line.find(element_search_key)
            part_keyword_flag = line.find(part_search_key)
            region_finisher_flag = line.find(region_finisher_key)
            end_keyword_flag = line.find(input_paragraph_end)
            if node_keyword_flag != -1:
                node_keyword_counter = node_keyword_counter + 1
                node_lines.append(line_counter)
            elif elem_keyword_flag != -1:
                elem_keyword_counter = elem_keyword_counter + 1
                elem_lines.append(line_counter)
            elif region_finisher_flag != -1:
                nset_keyword_counter = nset_keyword_counter + 1
                nset_lines.append(line_counter)
            elif part_keyword_flag != -1:
                part_keyword_counter = part_keyword_counter + 1
                part_lines.append(line_counter)
                file_identifiers.append(line.split(": ")[1].rstrip().replace(".", "_"))
            elif end_keyword_flag != -1:
                end_line_counter = line_counter
                break
            else:
                pass
        node_lines.append(end_line_counter)

    # Create a new dat file which will be used in the next steps by replacing commas with spaces and taking away the
    # LINE word at the beginning of some of the dat file lines

    new_earth_dat = fname[0:-4] + "_new.dat"
    with open(fname, "r+") as read_dat:
        with open(new_earth_dat, "w") as new_read_dat:
            for line in read_dat:
                replacements = line.replace(',', ' ').replace('LINE', ' ')
                new_read_dat.write(replacements)

    return node_lines, elem_lines, nset_lines, file_identifiers, new_earth_dat
