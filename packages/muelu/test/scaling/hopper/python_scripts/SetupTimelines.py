#!/usr/bin/env python
# Different setup times output
# Simply comment/uncomment the timers you are interested in
LABELS    = []
TIMELINES = []

# ========== TentativePFactory ==========
LABELS.append('tent P_1'); TIMELINES.append('TentativePFactory.*Build (level=1)')
# LABELS.append('tent P_2'); TIMELINES.append('TentativePFactory.*Build (level=2)')
# LABELS.append('tent P_3'); TIMELINES.append('TentativePFactory.*Build (level=3)')
# LABELS.append('tent P_4'); TIMELINES.append('TentativePFactory.*Build (level=4)')

# ========== SaPFactory ==========
LABELS.append('smoo P_1'); TIMELINES.append('SaPFactory.*Prolongator smoothing (level=1)')
# LABELS.append('smoo P_2'); TIMELINES.append('SaPFactory.*Prolongator smoothing (level=2)')
# LABELS.append('smoo P_3'); TIMELINES.append('SaPFactory.*Prolongator smoothing (level=3)')
# LABELS.append('smoo P_4'); TIMELINES.append('SaPFactory.*Prolongator smoothing (level=4)')

# ========== TransPFactory ==========
LABELS.append('R_1'); TIMELINES.append('TransPFactory.*Transpose P (level=1)')
# LABELS.append('R_2'); TIMELINES.append('TransPFactory.*Transpose P (level=2)')
# LABELS.append('R_3'); TIMELINES.append('TransPFactory.*Transpose P (level=3)')
# LABELS.append('R_4'); TIMELINES.append('TransPFactory.*Transpose P (level=4)')

# ========== RAPFactory ==========
LABELS.append('Ac_1'); TIMELINES.append('RAPFactory.*Computing Ac (level=1)')
# LABELS.append('Ac_2'); TIMELINES.append('RAPFactory.*Computing Ac (level=2)')
# LABELS.append('Ac_3'); TIMELINES.append('RAPFactory.*Computing Ac (level=3)')
# LABELS.append('Ac_4'); TIMELINES.append('RAPFactory.*Computing Ac (level=4)')

# ========== CoalesceDropFactory ==========
LABELS.append('CD_0'); TIMELINES.append('CoalesceDropFactory.*Build (level=0)')
# LABELS.append('CD_1'); TIMELINES.append('CoalesceDropFactory.*Build (level=1)')
# LABELS.append('CD_2'); TIMELINES.append('CoalesceDropFactory.*Build (level=2)')
# LABELS.append('CD_3'); TIMELINES.append('CoalesceDropFactory.*Build (level=3)')
# LABELS.append('CD_4'); TIMELINES.append('CoalesceDropFactory.*Build (level=4)')

# ========== ZoltanInteface ==========
# LABELS.append('Zoltan_1'); TIMELINES.append('Zoltan.Interface.*Build (level=1)')
LABELS.append('Zoltan_2'); TIMELINES.append('Zoltan.Interface.*Build (level=2)')
# LABELS.append('Zoltan_3'); TIMELINES.append('Zoltan.Interface.*Build (level=3)')
# LABELS.append('Zoltan_4'); TIMELINES.append('Zoltan.Interface.*Build (level=4)')
