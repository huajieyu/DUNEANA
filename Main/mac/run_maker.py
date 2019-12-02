import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()

maker.SetInputFile("/pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc0to30.root") # Run1 After Neutrino
maker.SetInputFile_add("/pnfs/dune/persistent/users/jiangl/pionana_1128/pionana_mc31to60.root") # Run1 After Neutrino
maker.SetOutputFile("output_test.root")

#maker.SetEntris(-1)
maker.MakeFile()

