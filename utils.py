# coding: utf-8
from __future__ import unicode_literals


def format_HTML(cols_formats):
    '''
    Columns for HTML table from DataFrame and their formats.
    For use in the to_html method of DataFrame.
    Parameter
        cols_formats: list of column names, formating pairs, e.g.,
        [['col1', '{:.6f}'], ['col2', '{:,.0f}' ]]
    Returns
        dictionary containing columns and lambdas
        that contain format strings, e.g.,
        {'columns': ['col1', 'col2'],
         'formatters: [lambda x: '{:.6f}'.format(x), lambda x: '{:,.0f}'.format(x)]
    '''
    columns = [e[0] for e in cols_formats]
    formatters = dict(zip(columns, [
                      eval(''.join(['lambda x: "', e[1],
                                    '".format(x)'])) for e in cols_formats]))
    return {'columns': columns,
            'formatters': formatters}
