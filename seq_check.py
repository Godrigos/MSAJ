#!/usr/bin/env python3

"""
This file is part of MSAJ.

MSAJ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MSAJ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MSAJ.  If not, see <https://www.gnu.org/licenses/>.

Copyright 2018 Rodrigo Aluizio
"""

from collections import defaultdict

def seq_check(msa):
    """
    Function to compare alignment files and check if all the strains
    names match between files.
    """
    id_list = defaultdict(list)
    i = 0

    for i in range(len(msa)):
        for j in range(len(msa[i])):
            id_list[i].append(msa[i][j].id)

    id_list = list(id_list.values())
    
    diff_ids = (set(id_list[0]).union(*id_list[1:]) -
                set(id_list[0]).intersection(*id_list[1:]))

    return diff_ids

