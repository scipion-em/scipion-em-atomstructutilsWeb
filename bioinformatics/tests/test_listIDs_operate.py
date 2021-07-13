# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_LISTIDs_OPERATE.PY
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from pathlib import Path
from pyworkflow.tests import *
from pyworkflow.protocol import *
from pwem.protocols.protocol_import import ProtImportPdb
from bioinformatics.protocols import ProtBioinformaticsListIDOperate as LOperate
from bioinformatics.protocols import ProtBioinformaticsDali as DALI


class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet("ligandLibraries")


class TestCreateOutput(TestImportBase):
    """ Create a realistic output of Dali using an old job (pdb id used: 5xjh)
    """

    def _importPDB(self, path):
        inputPdbData = 1 #file
        args = {'inputPdbData': inputPdbData,
                'pdbFile': path
                }

        prot1 = self.newProtocol(ProtImportPdb, **args)
        self.launchProtocol(prot1)
        global pdb
        pdb = prot1.outputPdb
        return pdb

    def getDalioutputs(self, protocol):
        fnBaseDir = self.dsModBuild.getFile(os.path.join("dali_output"))
        for fn in Path(fnBaseDir).rglob('*.txt'):
            protocol.constructOutput(str(fn), protocol)


class TestListOperate(TestCreateOutput):

    def test_1filter(self):
        """1. Filter a column using different ways
        """

        print("\n Filtering a SetOfDatabaseID object by a column")

        # Import a pdb file
        path = self.dsModBuild.getFile(os.path.join("proteinpdb", "5xjh.pdb"))
        pdb = self._importPDB(path)

        #Dali
        args = {'inputStructure': pdb,
                'method': 0, #search
                'title': 'test_dali_5xjh',
                'email': 'amparraperez@gmail.com'}

        prot1 = self.newProtocol(DALI, **args)
        prot1._store()
        self.getDalioutputs(prot1)
        prot1.setStatus(STATUS_FINISHED)

        global outputDali; global outputDali50; global outputDali25
        outputDali = prot1.outputDatabaseIds90
        outputDali50 = prot1.outputDatabaseIds50
        outputDali25 = prot1.outputDatabaseIds25

        self.assertIsNotNone(outputDali, "Error in creation of SetOfDatabaseID by DALI - It is NONE")
        self.assertTrue(outputDali.getSize() == 432, "The size of the SetOfDatabaseID is wrong. It should be 432 ")
        n_columns = len(list(outputDali.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 11, "The number of columns has changed. It should be 11")  # 2 fixed columns + 9 given


        # SET filtering :  _DaliZscore column ( >= 40)
        args = {'operation': 6,  # Filter
                'inputListID': outputDali,
                'filterColumn': '_DaliZscore',
                'filterOp': 2,  # >=
                'filterValue': 40,
                'removeDuplicates': True}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 7, "Error in filtering >= 40 of _DaliZscore column")


        # SET filtering :  _DaliZscore column ( 25 >= value <= 52 )
        args = {'operation': 6,  # Filter
                'inputListID': outputDali,
                'filterColumn': '_DaliZscore',
                'filterOp': 6,  # between
                'filterValue': 52,
                'filterValue2': 25,
                'removeDuplicates': True}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize()==8, "Error in filtering 25 >= value <= 52 of _DaliZscore column")

        # SET filtering :  _pdbId column ( startwith 6)
        args = {'operation': 6,  # Filter
                'inputListID': outputDali,
                'filterColumn': '_pdbId',
                'filterOp': 7,  # startwith
                'filterValue': "6",
                'removeDuplicates': True}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 57, "Error in filtering >_pdbId column. Startwith does not work")

        # SET filtering :  _DaliDescription column ( contains ESTERASE)
        args = {'operation': 6,  # Filter
                'inputListID': outputDali,
                'filterColumn': '_DaliDescription',
                'filterOp': 9,  # contains
                'filterValue': "ESTERASE",
                'removeDuplicates': True}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 104, "Error in filtering _DaliDescription column searching the ESTERASE word")


    def test_2Unique(self):
        """2. Keep unique entries regarding the DbId column
        """
        print("\n Keeping unique entries regarding the DbId column ")

        # Unique :  _pdbId column
        args = {'operation': 0,
                'inputListID': outputDali,
                'removeDuplicates': False
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432,
                        "Error in creation of a new SetOfDatabaseID. It has removed not unique entries")


    def test_3Union(self):
        """3. Union of different SetOfDatabaseID with remove duplicates or not
        """
        print("\n Union of 2 different SetOfDatabaseID with and without control of duplicates")

        args = {'operation': 1,
                'multipleInputListID': [outputDali50, outputDali25],
                'removeDuplicates': False
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 468, "Error in the union of 2 SetOfDatabaseID")


        # Removing duplicate entries

        args = {'operation': 1,
                'multipleInputListID':  [outputDali50, outputDali25],
                'removeDuplicates': True
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 320, "Error in the union of 2 SetOfDatabaseID")


    def test_4KeepColumns(self):
        """4. Keep only 2 columns and all entries
        """
        print("\n Keeping 2 interesting columns in a new SetOfDatabaseID")

        # Keep_column :  _pdbId column
        args = {'operation': 5,
                'inputListID': outputDali,
                'removeDuplicates': True,
                'keepColumns': '_pdbId ; _DaliZscore'
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432, "Error in creation of a new SetOfDatabaseID. It is different from the original regarding the number of entries")
        n_columns = len(list(setf.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 2+2, "Error. The number on columns is wrong. It should be 4") #2 fixed columns + 2 given


    def test_5Intersection(self):
        """6. Intersection between 2 different SetOfDatabaseID
        """
        print("\n Intersection between 2 different SetOfDatabaseID with control of duplicates")

        args = {'operation': 2,
                'inputListID': outputDali,
                'inputListID2': outputDali50,
                'removeDuplicates': True
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 311, "Error in the intersection of 2 SetOfDatabaseID")


    def test_6Difference(self):
        """6. Difference between 2 different SetOfDatabaseID
        """
        print("\n Difference between 2 different SetOfDatabaseIDs removing duplicates")

        args = {'operation': 3,
                'inputListID': outputDali,
                'inputListID2': outputDali50,
                'removeDuplicates': True
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 121, "Error in the difference of 2 SetOfDatabaseID")


    def test_7ChangeDbID(self):
        """7. Change the reference column in the SetOfDatabaseIDs
        """
        print("\n Change the reference column in the SetOfDatabaseIDs")

        args = {'operation': 4,
                'inputListID': outputDali,
                'newDb': "uniprot",
                'newDbId': "_DaliDescription",
                'removeDuplicates': True
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 302, "Error in the difference of 2 SetOfDatabaseID")