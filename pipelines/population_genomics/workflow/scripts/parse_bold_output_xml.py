# Print each BOLD match's species identification and similarity score, one per line.
# usage: parse_bold_output_xml.py XML_FILE

import sys
import xml.etree.ElementTree as ET

xml_file = sys.argv[1]
#out_tsv = sys.argv[2]

tree = ET.parse(xml_file)
root = tree.getroot()

for match in root:
	print(match.find("./taxonomicidentification").text.replace(" ", "_"), "\t", match.find("./similarity").text)

