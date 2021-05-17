#
# This file contains custom error classes
# This is intended to be used inside the parsers
# to provide additional information if something fails
#
# If additional more specific classes needed, subclass from ParserError
#


class ParserError(Exception):
    def __init__(self, parser_name, message):
        self.parser_name = parser_name
        self.message = message

    def __str__(self):
        return 'Error found while parsing output using "{}" parser: {}'.format(self.parser_name, self.message)
