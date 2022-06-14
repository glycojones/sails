import unittest
import os
import shutil
import sys
import sails
import test_data

class Test(unittest.TestCase):

    def setUp(self):
        '''
        Always run at start of test, e.g. for creating directory to store
        temporary test data producing during unit test
        '''
        self.test_data_path = os.path.dirname(test_data.__file__)
        assert os.path.exists(self.test_data_path)
        self.test_output = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'test_output')
        if not os.path.exists(self.test_output):
            os.makedirs ( self.test_output )


    def tearDown(self):
        '''
        Always run at end of test, e.g. to remove temporary data
        '''
        #shutil.rmtree(self.test_output)


    def test_database (self, verbose=True):
        '''
        Test Sails's fingerprint database
        '''

        print ("Testing Sails database")

        assert ( sails.number_of_fingerprints () == 0 )
        sails.initialise_fingerprints()
        print ("Now initialising fingerprints...")
        assert ( sails.number_of_fingerprints () > 0 )


if __name__ == '__main__':
    unittest.main()
