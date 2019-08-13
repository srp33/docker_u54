import glob
import sys
import os
import errno

test_dir = sys.argv[1]
write_access_required = sys.argv[2] == "True"
is_required = sys.argv[3] == "True"

if not is_required and not os.path.exists(test_dir):
    sys.exit(0)

try:
    for file_path in glob.glob(test_dir + "/*"):
        reading_file = open(file_path)
        reading_file.close()
except IOError as err:
    if err.errno == errno.EACCES:
        print(''.join(["ERROR: ", test_dir, " does not have read permissions"]))
    sys.exit(1)

if write_access_required:
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
