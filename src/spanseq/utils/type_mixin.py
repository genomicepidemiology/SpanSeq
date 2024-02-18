

class TypeValidation:

    @staticmethod
    def multiple_types(types, value, parser=None):
        print(value)
        for type in types:
            if not isinstance(value, type):
                if parser is None:
                    raise TypeError("""Only formats {} are acceted. {} was"""
                                    """ used""".format(", ".join(types), value)
                                    )
                else:
                    parser.error("""Only formats {} are acceted. {} was"""
                                    """ used""".format(", ".join(types), value)
                                    )
