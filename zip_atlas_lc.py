#!/usr/bin/env python
"""
@author: Sofia Rest
"""

import sys, argparse, configparser, zipfile, os

class zip_atlas_lc():
	def __init__(self):
		self.output_dir = None

	def define_args(self, parser=None, usage=None, conflict_handler='resolve'):
		if parser is None:
			parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)
		
		parser.add_argument('tnsnames', nargs='+', help='TNS names of the objects to download from ATLAS')
		parser.add_argument('--cfg_filename', default='atlas_lc_settings.ini', type=str, help='file name of ini file with settings for this class')

		return parser

	def set_output_dir(self, args):
		cfg = configparser.ConfigParser()
		try:
			print(f'Loading config file at {args.cfg_filename}')
			cfg.read(args.cfg_filename)
		except Exception as e:
			raise RuntimeError(f'ERROR: Could not load config file at {args.cfg_filename}!')

		self.output_dir = cfg['Input/output settings']['output_dir']
		print(f'Output directory: {self.output_dir}')

	def get_in_dirname(self, tnsname):
		return f'{self.output_dir}/{tnsname}'

	def get_out_filename(self, tnsname):
		return f'{self.output_dir}/SN{tnsname}'

	def getfilesfrom(self, directory):
		return filter(lambda x: not os.path.isdir(os.path.join(directory, x)), os.listdir(directory))

	def zip_directory(self, out_filename, in_dirname):
		zf = zipfile.ZipFile(f'{out_filename}.zip', mode='a', compression=zipfile.ZIP_DEFLATED)
		filestozip = self.getfilesfrom(in_dirname)
		for afile in filestozip:
			if '.DS_Store' in afile:
				continue
			print(f'# Adding {afile}')
			zf.write(os.path.join(in_dirname, afile), afile)
		zf.close()

	def zip_sne(self):
		args = self.define_args().parse_args()
		self.set_output_dir(args)

		for obj_index in range(0,len(args.tnsnames)):
			in_dirname = self.get_in_dirname(args.tnsnames[obj_index])
			print(f'Target directory: {in_dirname}')
			out_filename = self.get_out_filename(args.tnsnames[obj_index])
			print(f'Output file name: {out_filename}')
			self.zip_directory(out_filename, in_dirname)

if __name__ == "__main__":
	zip_atlas_lc = zip_atlas_lc()
	zip_atlas_lc.zip_sne()