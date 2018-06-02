import sys
'''
This script reads every line in the input file and merges
these in a string seperated by a specific delimiter given
as second command line argument. Additionally a third
(arbitrary) command line parameter can be given indictaing
whether the elements in the string should be enquoted.
'''
file_name = sys.argv[1]
seperator = sys.argv[2]

with open(file_name) as in_file:
    wanted = set()
    for line in in_file:
        line = line.strip()
        if line:
            if len(sys.argv) > 3:
                wanted.add("id:\"" + line + "\"")
            else:
                wanted.add(line)
print(seperator.join(wanted))
