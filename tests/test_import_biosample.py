import os
import sys

sys.path.insert(0, "../src/components")

import unittest
import pandas as pd
from importBioSample import ImportBioSample


class TestImportBioSample(unittest.TestCase):
    def setUp(self):
        self.file = "testbiosample.fastq"
        self.import_biosample = ImportBioSample(self.file)

    def test_importBioSample(self):
        resultDict, resultDf = self.import_biosample.importBioSample()
        fileExtension = os.path.splitext(self.file)[-1]
        self.assertIsInstance(resultDict, dict)
        self.assertIsInstance(resultDf, pd.DataFrame)
        self.assertTrue(fileExtension == ".fastq")
        if resultDict:
            first_key = list(resultDict.keys())[0]
            self.assertIn("sequence", resultDict[first_key])
            self.assertIn("quality", resultDict[first_key])
        if not resultDf.empty:
            self.assertIn("id", resultDf.columns)
            self.assertIn("sequence", resultDf.columns)
            self.assertIn("quality", resultDf.columns)


if __name__ == "__main__":
    unittest.main()
