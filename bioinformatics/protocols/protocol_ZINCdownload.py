# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
Describe your python module here:

"""
#IMPORT
import requests
import re
import warnings
import os
import glob

from pyworkflow.protocol import Protocol, params, Integer
import pyworkflow.object as pwobj
from pyworkflow.utils import Message

from bioinformatics.objects import SetOfSmallMolecules, SmallMolecule
from bioinformatics import Plugin



class ZINC_download(Protocol):
    """

    """
    _label = 'ZINC15 Download'
    _zinc_urlbase = 'http://zinc15.docking.org/'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label="General Download")

        form.addParam('zinc_mode', params.EnumParam, choices=['ZINC ID', 'List of ZINC IDs'],
                      label='Download a',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Download a single molecule from ZINC15 database using a proper ID (ZINC<number>) or using a'
                           ' list of these IDs. \n\n '
                           '* <number> is the number that identifies the compound (e.g. ZINC121 or ZINC000000000121 '
                           '(code of 12 digits)) ')
        form.addParam('zinc_id', params.StringParam, condition='zinc_mode == 0',
                      label='ZINC ID',
                      help='Download a single molecule from ZINC15 '
                           '(e.g. ZINC13108132 or ZINC000013108132 (hydroxytyrosol))')
        form.addParam('zinc_list', params.PathParam, condition='zinc_mode == 1',
                      label='ZINC List',
                      help='Download a list (.txt) of molecules from ZINC15 separates by comma, semicolon or in'
                           ' different lines.\n\n'
                           'Examples: \n\n '
                           '1- ZINC13108132,ZINC000013108133,ZINC000000000123 \n\n'
                           '2- ZINC13108132;ZINC000013108133;ZINC000000000123 \n\n'
                           '3- ZINC13108132\n   ZINC000013108133\n   ZINC000000000123 \n\n')
        form.addParam('format_general', params.EnumParam,
                      choices=['smi', 'sdf', 'mol2'], default=2,
                      label='Format',
                      help='You may download and create a set of small molecules from Smiles (.smi), '
                           'SDF (.sdf) and Tripos Mol2 (.mol2). The Smiled representation is in 2D.')


        form.addSection(label="Predefined subsets")
        form.addParam("weight_pol", params.BooleanParam,
                      label='Filter ZINC15 database using weight and polarity filters',
                      default=False, important=True,
                      help='Filter ZINC database of substances using different filters. \n\n'
                           'Warning! The size of molecules can be oversized')

        Tranches = form.addGroup("Weight and Polarity", condition='weight_pol')
        Tranches.addParam('weight', params.BooleanParam,
                          label='Weight', default=False,
                          important=True,
                          help='')

        form.addParam("completedb", params.BooleanParam,
                      label='Filter using others filters', default=False,
                      important=True,
                      help='Filter substances of ZINC database using different filters. \n\n'
                           'Warning! The size of molecules can be oversized')

        Availability = form.addGroup("Availability", condition='completedb')
        Availability.addParam('forSale', params.BooleanParam,
                              label='For Sale', default= False,
                              important=True,
                              help='Yes (for sale):  We download all the molecules that are for sale.'
                                   'It includes molecules in-stock, on-demand and boutique. \n\n'
                                   'No (not for sale):  We download all compounds that can '
                                   'not be bought as far as ZINC15 knows. ')
        Availability.addParam('available', params.EnumParam,
                              choices=['Options','Now', 'On-demand', 'Both'],
                              important=True,
                              label='Available', default=0,
                              condition='forSale',
                              help='Now means immediate delivery. It includes molecules in-stock and'
                                   ' have an agent.\n\n'
                                   'On demand: All substances that are for sale but not right now.\n\n'
                                   'Both: Compounds you can get in 8-10 weeks at modest prices. '
                                   'It includes molecules in-stock, agent and on-demand. \n\n'
                                   'Options: No filter (Now or On demand) will be applied')
        Availability.addParam('bb', params.BooleanParam,
                              label='Preparative Quantities available', default=False,
                              condition='forSale',
                              help='Available in preparative quantities, typically at least 250 mg.')
        Availability.addParam('boutique', params.BooleanParam,
                              label='Expensive substances', default=False,
                              condition='forSale',
                              help='Boutique substances are generally much more expensive than $100/sample, '
                                   'often made to order, yet still cheaper than making it yourself.')


        Bioactive = form.addGroup("Bioactive in and approved by", condition='completedb')
        Bioactive.addParam('fda', params.EnumParam,
                           choices=['Options','FDA', 'World but not FDA', 'Major jurisdictions', "Investigational"],
                           important=True,
                           label='Approved by', default=0,
                           help='FDA: FDA Approved drugs, per DrugBank.\n\n'
                                'World but not FDA: Drugs approved, but not by the FDA.\n\n'
                                'Major jurisdictions: Approved drugs in major jurisdictions, including the FDA,'
                                ' i.e DrugBank approved.\n\n'
                                'Investigational: Investigational compounds - in clinical trials - not approved '
                                'or used as drugs.\n\n'
                                'Options: No filter from the above will be applied ')
        Bioactive.addParam('vitro', params.BooleanParam,
                           label='In vitro', default=False,
                           help='Substances reported or inferred active at 10 uM or better in direct binding assays')
        Bioactive.addParam('cells', params.BooleanParam,
                           label='In Cells', default=False,
                           help='Substances reported or inferred active in cells')
        Bioactive.addParam('vivo', params.BooleanParam,
                           label='In vivo', default=False,
                           help='Substances tested in animals including man')
        Bioactive.addParam('notman', params.BooleanParam,
                           label='Exclude man', default=False,
                           condition= 'vivo',
                           help='Substances tested in animals but not in man, e.g. DrugBank Experimental')
        Bioactive.addParam('man', params.BooleanParam,
                           label='In Human', default=False,
                           help='Substances that have been in man')
        Bioactive.addParam('human', params.BooleanParam,
                           label='In Human but not approved', default=False,
                           condition='man',
                           help='Substances that have been in man, but not approved or in trials, e.g'
                                ' nutriceuticals and many metabolites')

        Biogenic = form.addGroup("Metabolites and Natural products", condition='completedb')
        Biogenic.addParam('biogenic', params.BooleanParam,
                          label='Biogenic', default=False,
                          important=True,
                          help='Made by nature, including primary metabolites (metabolites) and secondary '
                               'metabolites (natural products)')
        Biogenic.addParam('metabolites', params.BooleanParam,
                          label='Metabolites', default=False,
                          condition='biogenic',
                          help='Primary metabolites of any species, including humans')
        Biogenic.addParam('type_metabolites', params.EnumParam,
                           choices=['Human', 'Not Human ', 'Both'],
                           label='Type of 1º metabolites', default=2,
                           condition='biogenic and metabolites',
                           help='Human: metabolites observed in humans.\n\n'
                                'Not Human: Primary metabolites - also known as metabolites, not reported in '
                                'humans.\n\n'
                                'Both: Primary metabolites - also known as metabolites, not reported in humans. ')
        Biogenic.addParam('natural_products', params.BooleanParam,
                          label='Natural products', default=False,
                          condition='biogenic',
                          help='Natural products, also known as secondary metabolites, i.e. '
                               'explicitly excluding metabolites')

        Reactivity = form.addGroup("Reactivity", condition='completedb')
        Reactivity.addParam('anodyne', params.BooleanParam,
                            label='Anodyne', default=False,
                            important=True,
                            help='Substances matching no reactivity patterns, including PAINS (pan-assay interference)')
        Reactivity.addParam('reactive', params.EnumParam,
                            choices=['All', 'Clean', 'Mild reactivity', 'Reactive', 'Hot or unstable'],
                            label='Reactive', default=0,
                            important=True,
                            condition='not anodyne',
                            help='Clean: Substances with "clean" reactivity, i.e. many PAINS (pan-assay interference)'
                                 ' patterns allowed.\n\n'
                                 'Mild or standard reactivity: Substances that are weakly reactive,'
                                 ' typically as a nucleophile or electrophile.\n\n'
                                 'Reactive: Substances matching one of the reactive patterns. \n\n'
                                 'Hot, unstable or irrelevant for screening: Substances matching one of the hot '
                                 'patterns. \n\n'
                                 'All: Reactive molecules in general, including PAINS, hot and ZINC12 "clean" filters.')

        Others = form.addGroup("Others", condition='completedb')
        Others.addParam('named', params.BooleanParam,
                         label='Named', default=False,
                         help='Compounds that have names in ZINC')
        Others.addParam('aggregators', params.BooleanParam,
                         label='Colloidal aggregates', default=False,
                         help='These compounds have been observed to form colloidal aggregates under assay conditions')

        form.addSection(label="Catalogs")

        form.addSection(label="Search similar compounds")

        form.addSection(label="ZINC12. Advanced Filters")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('generalDownload')
        self._insertFunctionStep('createOutputStep')


    def createOutputStep(self):

        outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(), suffix='_fromZINC')

        for file in glob.glob(self._getExtraPath("*")):

            smallMolecule = SmallMolecule(smallMolFilename=file)

            # Print the image of the small molecule
            if file.endswith(".pdb"):
                fnRoot = os.path.splitext(os.path.split(file)[1])[0]
                fnOut = self._getExtraPath("%s.png" % fnRoot)
                args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (file, fnOut)
                try:
                    Plugin.runRDKit(self, "python3", args)
                    smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                except:
                    smallMolecule._PDBLigandImage = pwobj.String("Not available")

            else:
                smallMolecule._PDBLigandImage = pwobj.String("Not available")

            outputSmallMolecules.append(smallMolecule)

        self._defineOutputs(SetOfSmallMolecules_ZINC=outputSmallMolecules)



    def generalDownload(self):
        """
        """

        # single ZINC ID
        if self.zinc_mode.get() == 0:
            zinc_id = self.zinc_id.get()

            if zinc_id =="": # Checks
                error = " ZINC ID is empty. Fill in the space designated for it"
                raise Exception(error)
            elif not re.search(r'^zinc?\d+$', zinc_id, re.IGNORECASE):
                error = " ZINC ID has a wrong format. Please check it "
                raise Exception(error)

            IDs = [zinc_id.replace(" ", "").upper()]


        # multiple ZINC ID. LIST
        elif self.zinc_mode.get() == 1:
            with open(self.zinc_list.get(), 'r') as f:
                file = f.read()
            if file.replace(" ", "").upper()[-1] == '\n':
                IDs = re.split(",|;|\n", file.replace(" ", "").upper()[:-1])
            else:
                IDs = re.split(",|;|\n", file.replace(" ", "").upper())

            if len(IDs)==0:
                error = "The list of ZINC IDs is empty"
                raise Exception(error)

        # Format
        if self.format_general.get() == 0:
            format_d = ".smi"
        elif self.format_general.get() == 1:
            format_d = ".sdf"
        else:
            format_d = ".mol2"

        discarded_ids = []
        for zinc_id in IDs:
            url = self._zinc_urlbase + 'substances/'
            url += zinc_id + format_d
            req = requests.get(url)


            if req.status_code != 200 and len(IDs) == 1:
                error =   "There was a request error in the download:\n" \
                          "    ZINC ID: %s" \
                          "    Status code: %s \n" \
                          "    Error reason: %s \n" \
                          "    URL: %s \n" \
                          "Check if ZINC ID is correct, please. " % (zinc_id, req.status_code, req.reason, req.url)
                raise Exception(error)

            elif req.status_code != 200:
                warning = "There was a request error in the download:\n" \
                          "    ZINC ID: %s \n" \
                          "    Status code: %s \n" \
                          "    Error reason: %s \n" \
                          "    URL: %s \n" \
                          "Check if ZINC ID is correct, please. " \
                          "The ZINC ID %s has been ignored.\n" % (zinc_id, req.status_code, req.reason, req.url, zinc_id)
                discarded_ids.append(zinc_id)
                warnings.warn(warning) #Warming
            else:
                file = self._getExtraPath(self.getFilename_fromHeader(req.headers)) #Create molecule file in extra folder
                with open(file, "wb") as f:
                    f.write(req.content)


        # Summary
        print('\nNUMBER OF DISCARDED FILES -->  ', len(discarded_ids))
        print('\n'.join('{}: {}'.format(k+1, discarded_ids[k]) for k in range(len(discarded_ids))))
        print('\nNUMBER OF TOTAL DOWNLOADED FILES -->  ', len(os.listdir(self._getExtraPath()+'/')))
        list_dirfiles = os.listdir(self._getExtraPath() + '/')
        print('\n'.join('{}: {}'.format(k+1, list_dirfiles[k]) for k in range(len(list_dirfiles))))
        print('\n')


    # --------------------------- Utils functions -----------------------------------
    def getFilename_fromHeader(self, header):
        """
        Get filename from header dictionary and content-disposition key
        """
        content = header.get('content-disposition')
        if not content:
            return None
        filename = re.findall('filename=(.+)', content)
        if len(filename) == 0:
            return None
        return filename[0]


    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        #if self.isFinished():
        #    summary.append("This protocol has printed *%s* %i times." % (self.message, self.times))
        #return summary
        pass

    def _methods(self):
        #methods = []

        #if self.isFinished():
        #    methods.append("%s has been printed in this run %i times." % (self.message, self.times))
        #    if self.previousCount.hasPointer():
        #        methods.append("Accumulated count from previous runs were %i."
        #                       " In total, %s messages has been printed."
        #                       % (self.previousCount, self.count))
        #return methods
        pass