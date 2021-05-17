from qcparsers.parsers import parser_optimization, parser_irc, parser_rasci
from qcparsers.parsers import parser_frequencies, parser_basic, parser_cis, parser_fchk
import unittest
import pickle


store = False


class ParsersTest(unittest.TestCase):

    def atest_optimization_1(self):

        with open('optimization_1.out', 'r') as f:
            qchem_output = f.read()

        data = parser_optimization(qchem_output)

        if store:
            with open('optimization_1.pkl', 'wb') as outfile:
                pickle.dump(data, outfile, protocol=2)

        with open('optimization_1.pkl', 'rb') as stream:
            data_ref = pickle.load(stream)

        self.assertDictEqual(data, data_ref)

    def atest_rasci_1(self):

        with open('rasci_1.out', 'r') as f:
            qchem_output = f.read()

        data = parser_rasci(qchem_output)

        if store:
            with open('rasci_1.pkl', 'wb') as outfile:
                pickle.dump(data, outfile, protocol=2)

        with open('rasci_1.pkl', 'rb') as stream:
            data_ref = pickle.load(stream)

        print(data)
        self.assertDictEqual(data, data_ref)


def add_test(cls, f_name, parser):
    def test_method(self):

        with open('{}.out'.format(f_name), 'r') as f:
            qchem_output = f.read()

        data = parser(qchem_output)

        if store:
            with open('{}.pkl'.format(f_name), 'wb') as outfile:
                pickle.dump(data, outfile, protocol=2)

        with open('{}.pkl'.format(f_name), 'rb') as stream:
            data_ref = pickle.load(stream)

        self.assertDictEqual(data, data_ref)

    # test_method.__doc__ = 'Test for file: {}.out'.format(t_name)
    test_method.__name__ = 'test_{}'.format(f_name)
    setattr(cls,test_method.__name__,test_method)


for f_name, parser in zip(['optimization_1', 'rasci_1', 'rasci_2', 'irc_1',
                           'frequencies_1', 'simple_1', 'cis_1', 'cis_2', 'fchk_1'],
                          [parser_optimization, parser_rasci, parser_rasci, parser_irc,
                           parser_frequencies, parser_basic, parser_cis, parser_cis, parser_fchk]):
    add_test(ParsersTest, f_name, parser)

