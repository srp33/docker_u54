import os
import sys
import textwrap
from yaml import load
from yaml import Loader

yaml_file_path = sys.argv[1]
tag = sys.argv[2]
out_file_path = sys.argv[3]
help_file_path = sys.argv[4]

def get_yaml_values(yaml_dict, key):
    for info in yaml_dict[key]:
        key = list(info.keys())[0]
        value = info[key]

        yield key, value

command = os.path.basename(yaml_file_path).replace(".yaml", "")

with open(yaml_file_path, 'r') as yaml_file:
    yaml_text = yaml_file.read()
yaml_dict = load(yaml_text, Loader=Loader)

output = "#! /bin/bash\n\n"

# Add comment header.
output += "#" * 60
output += "\n# This script was generated by {}.\n".format(sys.argv[0])
output += "# Changes should be made to yaml files.\n".format(sys.argv[0])
output += "#" * 60
output += "\n\n"

# Specify script settings and load functions.
output += "set -o errexit\n"
output += "source /jicama/helper/check_functions\n\n"

# Generate help output dynamically
help_header_template = "{}\n{}{}\n\n"
help_output = help_header_template.format("NAME", "\t", command)
help_output += help_header_template.format("DESCRIPTION", "\t", "\n\t".join(textwrap.wrap(yaml_dict["description"])))

help_options = ""
for arg, meta in get_yaml_values(yaml_dict, "args"):
    opts = meta["opts"].replace(" | ", ", ")

    parenthetical_statements = []
    if "default" in meta:
        parenthetical_statements.append("optional")
    if "allow_multiple" in meta and meta["allow_multiple"] == True:
        statement = "may be specified multiple times"
        if "min_required" in meta:
            statement += " [minimum: {}]".format(meta["min_required"])
        parenthetical_statements.append(statement)

    parenthetical_statement = ""
    if len(parenthetical_statements) > 0:
        parenthetical_statement = " (" + ", ".join(parenthetical_statements) + ")"

    description = meta["description"]
    if "example" in meta:
        description += " Example value: {}.".format(meta["example"])
    description = "\n\t\t".join(textwrap.wrap(description))

    help_options += "\t{}{}\n\t\t{}\n\n".format(opts, parenthetical_statement, description)
help_output += help_header_template.format("OPTIONS", "", help_options.rstrip())

help_volumes = "\t" + "\n\t".join(textwrap.wrap("When executing software within a Docker container, you can use volumes to map directories outside the container to directories inside the container. The -v argument is used for this. Multiple -v arguments can be specified when a container is executed. The value for each -v argument has two parts, separated by a colon. The first part should be an absolute path to a directory on the local (host) computer. The second part should be an absolute path to a directory within the container. To avoid problems with permissions, you should ensure that the directories on the host computer have been created before executing the container. This command requires the following volumes:")) + "\n\n"
for volume, meta in get_yaml_values(yaml_dict, "volumes"):
    help_volumes += "\t<absolute path on local computer>:/volumes/{}\n".format(volume)
    help_volumes += "\t\t" + "\n\t\t".join(textwrap.wrap(meta["description"])) + "\n"

    if "write_access" in meta and meta["write_access"] == True:
        help_volumes += "\t\tThe current user must have read/write permission on this directory.\n\n"
    else:
        help_volumes += "\t\tThe current user must have read permission on this directory.\n\n"

    if "required" in meta and meta["required"] == False:
        help_volumes += "\t\tThis volume is optional.\n\n"
help_output += help_header_template.format("VOLUMES", "", help_volumes.rstrip())

help_example = "\tdocker run \\\n"
count = 0
for volume, meta in get_yaml_values(yaml_dict, "volumes"):
    count += 1
    help_example += "\t  -v /my/local/path{}:/volumes/{} \\\n".format(count, volume)
help_example += "\t  --user \$(id -u):\$(id -g) \\\n"
help_example += "\t  --rm \\\n"
help_example += "\t  srp33/jicama:{} \\\n".format(tag)
help_example += "\t  {} \\\n".format(command)
for arg, meta in get_yaml_values(yaml_dict, "args"):
    if "example" in meta and "default" not in meta:
        example_opt = meta["opts"].split(" | ")[0]

        if meta["example"].startswith("\""):
            example_values = [meta["example"]]
        else:
            example_values = meta["example"].split(" ")

        example_output = example_opt + " " + " \\\n        {} ".format(example_opt).join(example_values)
        help_example += "\t    {} \\\n".format(example_output)
help_output += help_header_template.format("EXAMPLE", "", help_example.rstrip().rstrip("\\"))

# Save help file
with open(help_file_path, 'w') as help_file:
    help_file.write(help_output)

# Print help function
output += "function show_help {\n"
output += "  cat /jicama/docs/{}\n".format(os.path.basename(help_file_path))
output += "}\n\n"

# Declare a variable for each argument with its default value.
for arg, meta in get_yaml_values(yaml_dict, "args"):
    if "default" in meta:
        default = "\"" + str(meta["default"]) + "\""
    else:
        if "allow_multiple" in meta and meta["allow_multiple"] == True:
            default = "()"
        else:
            default = "\"Null\""

    output += "{}={}\n".format(arg, default)

output += "\n"

# Parse value from user for each argument.
output += "ARGNUM=$#\n\n"
output += "for (( i=1; i<=ARGNUM; i++ )); do\n"
output += "  OPTARG=$((i+1))\n"
output += "  case ${!i} in\n"

for arg, meta in get_yaml_values(yaml_dict, "args"):
    output += "    " + meta["opts"] + " )\n"
    output += "      check_args \"${!OPTARG}\" \"${!i}\" || exit 1\n"

    if "allow_multiple" in meta and meta["allow_multiple"] == True:
        output += "      " + arg + "+=(\"${!OPTARG}\")\n"
    else:
        output += "      " + arg + "=\"${!OPTARG}\"\n"

    output += "      i=$((i+1))\n"
    output += "      ;;\n"

# Print help behavior.
output += "    -h | --help )\n"
output += "      show_help\n"
output += "      exit 0\n"
output += "      ;;\n"
output += "    * )\n"
output += "      echo \"Invalid option: ${!i}\" 1>&2\n"
output += "      exit 1\n"
output += "      ;;\n"
output += "  esac\n"
output += "done\n\n"

# Check whether required arguments are present.
# An argument is required if there is no default value.
for arg, meta in get_yaml_values(yaml_dict, "args"):
    if "default" not in meta:
        if "allow_multiple" in meta and meta["allow_multiple"] == True:
            min_required = "1"
            if "min_required" in meta:
                min_required = str(meta["min_required"])

            output += "if [[ ${#" + arg + "[@]} -lt " + min_required + " ]]\n"
        else:
            output += "if [[ \"${" + arg + "}\" == \"Null\" ]]\n"

        output += "then\n"
        output += "  echo \"ERROR: The " + meta["opts"] + " argument must be provided.\"\n"
        output += "  echo \n"
        output += "  show_help\n"
        output += "  exit 1\n"
        output += "fi\n\n"

# Check whether argument values are acceptable.
for arg, meta in get_yaml_values(yaml_dict, "args"):
    if "acceptable_values" in meta:
        output += "if [[ "
        for value in meta["acceptable_values"]:
            output += "\"${" + arg + "}\" != \"" + str(value) + "\" && "
        output = output.rstrip("&& ")
        output += " ]]\n"
        output += "then\n"
        output += "  echo \"ERROR: The " + arg + " argument must be " + " or ".join([str(x) for x in meta["acceptable_values"]]) + ".\"\n"
        output += "  show_help\n"
        output += "  exit 1\n"
        output += "fi\n\n"

# Check for the necessary directories which are only created by volumes
output += "MISSING_VOLUMES=()\n"
output += "EXIT_CODE=()\n\n"

# Check for volumes
for volume, meta in get_yaml_values(yaml_dict, "volumes"):
    if "required" not in meta or meta["required"] == True:
        output += "[[ -d " + volume + " ]] || { MISSING_VOLUMES+=(" + volume + ") && EXIT_CODE=1; }\n"
output += "\n"

output += "if [[ ${EXIT_CODE} = 1 ]]\n"
output += "then\n"
output += "  echo \"The following volumes are missing: ${MISSING_VOLUMES[@]}\"\n"
output += "  show_help\n"
output += "  exit 1\n"
output += "fi\n\n"

# Check permissions of each directory
for volume, meta in get_yaml_values(yaml_dict, "volumes"):
    write_access = meta.get("write_access", False)
    required = meta.get("required", True)

    output += "python /jicama/helper/check_permissions.py " + volume + " " + str(write_access) + " " + str(required) + " || exit 1\n"

output += "\n"

#output += "set -o xtrace\n\n"

# Add the commands for this specific task
output += yaml_dict["command_template"] + "\n\n"

# Create a receipt with version numbers, etc.
output += "echo " + ("=" * 50) + "\n"
output += "echo Receipt\n"
output += "echo " + ("=" * 50) + "\n"
output += "echo\n"
output += "echo Timestamp: $(date '+%d/%m/%Y %H:%M:%S')\n"
output += "echo\n\n"
output += yaml_dict["receipt_commands"]

print("Saving to {}.".format(out_file_path))
with open(out_file_path, 'w') as out_file:
    out_file.write(output)
