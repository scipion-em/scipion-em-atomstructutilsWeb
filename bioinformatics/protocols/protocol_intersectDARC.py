
import os, re
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam, IntParam, FloatParam
from pyworkflow.object import Float, Integer
from bioinformatics.objects import SmallMolecule, SetOfSmallMolecules
from rosetta.objects import SetScores


class ProtBioinformaticsDARCPath(EMProtocol):
    """Retrieve information about paths of molecules with better score"""
    _label = 'Retrieve Path of Molecules'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set to filter:', allowsNull=False)
        form.addParam('secondSet', PointerParam, pointerClass="SetScores",
                       label='Set with Scores:', allowsNull=False)
        """
        form.addParam('filterColumn', StringParam,
                       label='Filter column:',
                       help='It must exist in the input object.')
        """


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('retrivesmallMolecules')

    def retrivesmallMolecules(self):

        outputSet = self.inputSet.get().create(self._getPath())

        # Get paths
        IDs = []

        for entry in self.secondSet.get():
            IDs.append((entry.getzincID()).upper())

        for entry in self.inputSet.get():
            original_path = entry.getFileName()
            basename = os.path.basename(original_path)
            ID = re.split("[_.]", basename)[0]

            if ID.upper() in IDs:
                outputSet.append(entry)


        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)