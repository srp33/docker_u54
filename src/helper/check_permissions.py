import glob
import sys
import os
import errno

test_dir = sys.argv[1]
needed_permissions = sys.argv[2]
#test_file = ""

if needed_permissions == "Read":
    try:
#        test_file = sys.argv[3]
#        reading_file = open('/'.join([test_dir, test_file]), 'r')
#        reading_file.close()
        for file_path in glob.glob(test_dir + "/*"):
            reading_file = open(file_path)
            reading_file.close()
    except IOError as err:
        if err.errno == errno.EACCES:
            print(''.join(["ERROR: ", test_dir, " does not have read permissions"]))
#        elif err.errno == errno.ENOENT:
#            print(''.join([test_file, " does not exist"]))
        sys.exit(1)
else:
    try:
        writing_file = open('/'.join([test_dir, "temp.txt"]), 'w')
        writing_file.write("TEST TEXT")
        writing_file.close()
    except IOError:
        print(''.join(["ERROR: ", test_dir, " does not have write permissions"]))
        sys.exit(1)
    try:
        reading_file = open('/'.join([test_dir, "temp.txt"]), 'r')
        reading_file.close()
        os.remove('/'.join([test_dir, "temp.txt"]))
        sys.exit(0)
    except IOError:
        print(''.join(["ERROR: ", test_dir, " does not have read permissions"]))
        sys.exit(1)
