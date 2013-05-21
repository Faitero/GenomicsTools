#!/usr/bin/python 
import subprocess;
import logging;
import glob;
import re;

class Sample:
	
	
	def __init__(self, name):
		self.name = name
		self.attribute = {}

	def setAttrib(self, name, value):
		self.attribute[name] = value

	def getAttrib(self, name):
		return self.attribute[name]



def parse_bam_stat_tab(samples, pattern="*bam_stat.tab" ):
	print pattern
	bam_stat_tab = {}
	for sample in samples:
		bam_stat_tab[sample] = glob.glob(sample+"/"+pattern)
		logging.debug(bam_stat_tab[sample])
		with open(bam_stat_tab[sample][0], 'r') as f:
			
			for line in f:
				print line
				# match = re.match(r'^(.+)[\t](\d+)', line)

				match = re.match(r'(.+)\s(\d+)$', line)
				if match:
					logging.debug(line+"Matched!")
					logging.debug( match.group(1) )
					re.sub(r' ','_', match.group(1))
					logging.debug( match.group(1) )
					logging.debug( match.group(2) )

				else:
					print("No Match")




def get_samples(pattern="Tophat*"):
	logging.debug(pattern);
	# tab_files = subprocess.check_output(['ls', pattern]).splitlines();
	tab_files = glob.glob(pattern);
	samples = {}
	for sample in tab_files:
		samples[sample] = Sample(sample)
	return samples


def main():
	logging.getLogger(__name__)
	# logging.basicConfig(filename='buildQCTable.log', level=logging.DEBUG)
	logging.basicConfig( level=logging.DEBUG)
	logging.basicConfig(format='%(asctime)s %(message)s')
	logging.info("Getting Samples")
	samples = get_samples()
	for key in samples:
		logging.debug( key)

	parse_bam_stat_tab(samples)


if __name__ == '__main__':
	main()