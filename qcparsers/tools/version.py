#
# Here the functions related to Q-Chem version handling system
#

def get_version_output(output):
    """
    Obtain the version from a Q-Chem output

    :param output: Q-Chem ouput in plain text
    :return: the version
    """

    class QChemVersion:
        def __init__(self, string):

            string_version = string.split()[1]
            self._major = string_version.split('.')[0]
            self._minor = string_version.split('.')[1]

            string_branch = string.split()[2]
            self._devel = True if '(devel)' in string_branch else False

        def __str__(self):
            dev = 'dev' if self.is_development else ''
            return '{}.{} {}'.format(self.major, self.minor, dev)

        def __eq__(self, other):

            """
            Here put the logic for more sophisticated comparison
            between versions

            :param other: string/QchemVersion
            :return:
            """

            if isinstance(other, QChemVersion):
                return True if self.__str__() == other.__str__() else False

            o_major = other.split('.')[0]
            o_minor = other.split('.')[1]

            if int(o_major) == self.major:

                # handle expresions like 2.3+
                if '+' in o_minor[-1]:
                    if self.minor >= int(o_minor[:-1]):
                        return True
                    else:
                        return False

                if int(o_minor) == self.minor:
                    return True

            return False

        @property
        def major(self):
            return int(self._major)

        @property
        def minor(self):
            return int(self._minor)

        @property
        def is_development(self):
            return self._devel

    index = output[:500].find('\n Q-Chem')
    string = output[index: index + 30]

    return QChemVersion(string)


def get_compatibility_list_from_parser(parser):
    """
    Extract the docstring information (if any) about compatible versions

    :param parser: a parser function
    :return: list of versions
    """
    docstring = parser.__doc__
    if docstring is None:
        return None

    lines = docstring.split('\n')
    for line in lines:
        if 'compatibility' in line.lower():
            try:
                return [version.strip() for version in line.split(':')[1].split(',')]
            except IndexError:
                continue

    return None


if __name__ == '__main__':

    from qcparsers.parsers import parser_basic

    with open('../../test/cis_1.out', 'r') as f:
        output = f.read()

    version = get_version_output(output)
    compatibility_list = get_compatibility_list_from_parser(parser_basic)

    if version in compatibility_list:
        print('Parser compatible')
    else:
        print('Parser not compatible')
