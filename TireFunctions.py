import numpy as np

L = dict(
    LFZO=1.0,
    LCZ=1.0,
    LCX=1.0, LMUX=1.0, LEX=1.0, LKX=1.0, LHX=1.0, LVX=1.0, LGAX=1.0,
    LCY=1.0, LMUY=1.0, LEY=1.0, LKY=1.0, LHY=1.0, LVY=1.0, LGAY=1.0,
    LGAZ=1.0, LTR=1.0, LRES=1.0, LVMX=1.0, LMY=1.0,
    LXAL=1.0, LYKA=1.0, LVYKA=1.0, LS=1.0
)

FY_Params = ["PCY1", "PDY1", "PDY2", "PDY3", "PEY1", "PEY2","PEY3","PEY4",
             "PKY1", "PKY2", "PKY3", "PHY1", "PHY2", "PHY3",
             "PVY1", "PVY2", "PVY3", "PVY4"]

MZ_Params = ["QBZ1", "QBZ2", "QBZ3", "QBZ4", "QBZ5", "QBZ9", "QBZ10", "QCZ1",
             "QDZ1", "QDZ2", "QDZ3", "QDZ4", "QDZ6", "QDZ7", "QDZ8", "QDZ9",
             "QEZ1", "QEZ2", "QEZ3", "QEZ4", "QEZ5", "QHZ1", "QHZ2", "QHZ3", "QHZ4"]

FX_Params = ["PCX1", "PDX1", "PDX2", "PDX3", "PEX1", "PEX2", "PEX3", "PEX4",
             "PKX1", "PKX2", "PKX3", "PHX1", "PHX2", "PVX1", "PVX2"]

def FY(X, p, fz0):
    SA, FZ, IA = X
    ALPHA = np.radians(SA)
    GAMMA = np.radians(IA)

    (PCY1, PDY1, PDY2, PDY3, PEY1, PEY2,PEY3,PEY4,
     PKY1, PKY2, PKY3, PHY1, PHY2, PHY3,
     PVY1, PVY2, PVY3, PVY4) = p

    FZ0 = fz0 * L["LFZO"]
    DFZ = (FZ - FZ0) / FZ0
    GAMY = GAMMA * L["LGAY"]

    SHY = (PHY1 + PHY2 * DFZ) * L["LHY"] + PHY3 * GAMY
    ALPHAY = ALPHA + SHY
    SVY = FZ * ((PVY1 + PVY2 * DFZ) * L["LVY"] + (PVY3 + PVY4 * DFZ) * GAMY) * L["LMUY"]

    CY = PCY1 * L["LCY"]
    MUY = (PDY1 + PDY2 * DFZ) * (1.0 - PDY3 * GAMY ** 2) * L["LMUY"]
    DY = MUY * FZ
    KY = (PKY1 * fz0 * np.sin(2.0 * np.arctan(FZ / (PKY2 * fz0 * L["LFZO"]))) * (1.0 - PKY3 * np.abs(GAMY)) * L['LFZO'] * L['LKY'])
    BY = KY / (CY * DY + 1e-9)
    EY = (PEY1 + PEY2 * DFZ) * (1.0 - (PEY3 + PEY4 * GAMY) * np.sign(ALPHAY)) * L['LEY']

    return DY * np.sin(CY * np.arctan(BY * ALPHAY - EY * (BY * ALPHAY - np.arctan(BY * ALPHAY)))) + SVY

def MZ(X, q, p, fz0, r0):
    SA, FZ, IA = X
    ALPHA = np.radians(SA)
    GAMMA = np.radians(IA)

    (PCY1, PDY1, PDY2, PDY3, PEY1, PEY2, PEY3, PEY4,
     PKY1, PKY2, PKY3, PHY1, PHY2, PHY3,
     PVY1, PVY2, PVY3, PVY4) = p

    (QBZ1, QBZ2, QBZ3, QBZ4, QBZ5, QBZ9, QBZ10, QCZ1,
     QDZ1, QDZ2, QDZ3, QDZ4, QDZ6, QDZ7, QDZ8, QDZ9,
     QEZ1, QEZ2, QEZ3, QEZ4, QEZ5, QHZ1, QHZ2, QHZ3, QHZ4) = q

    FZ0 = fz0 * L['LFZO']
    DFZ = (FZ - FZ0) / FZ0
    GAMY = GAMMA * L['LGAY']
    GAMZ = GAMMA * L['LGAZ']
    SHY = (PHY1 + PHY2 * DFZ) * L['LHY'] + PHY3 * GAMY
    SVY = FZ * ((PVY1 + PVY2 * DFZ) * L['LVY'] + (PVY3 + PVY4 * DFZ) * GAMY) * L['LMUY']
    ALPHAY = ALPHA + SHY
    KY = (PKY1 * fz0 * np.sin(2.0 * np.arctan(FZ / (PKY2 * fz0 * L['LFZO']))) * (1.0 - PKY3 * np.abs(GAMY)) * L['LFZO'] * L['LKY'])
    CY = PCY1 * L['LCY']
    MUY = (PDY1 + PDY2 * DFZ) * (1.0 - PDY3 * GAMY ** 2) * L['LMUY']
    DY = MUY * FZ
    BY = KY / (CY * DY + 1e-9)
    EY = (PEY1 + PEY2 * DFZ) * (1.0 - (PEY3 + PEY4 * GAMY) * np.sign(ALPHAY)) * L['LEY']
    FY0 = DY * np.sin(CY * np.arctan(BY * ALPHAY - EY * (BY * ALPHAY - np.arctan(BY * ALPHAY)))) + SVY

    SHT = QHZ1 + QHZ2 * DFZ + (QHZ3 + QHZ4 * DFZ) * GAMZ
    ALPHAT = ALPHA + SHT
    SHF = SHY + SVY / (KY + 1e-9)
    ALPHAR = ALPHA + SHF

    BT = ((QBZ1 + QBZ2 * DFZ + QBZ3 * DFZ ** 2) * (1.0 + QBZ4 * GAMZ + QBZ5 * np.abs(GAMZ)) * L['LKY'] / L['LMUY'])
    CT = QCZ1
    DT = (FZ * (QDZ1 + QDZ2 * DFZ) * (1.0 + QDZ3 * GAMZ + QDZ4 * GAMZ ** 2) * (r0 / fz0) * L['LTR'])
    ET = ((QEZ1 + QEZ2 * DFZ + QEZ3 * DFZ ** 2) * (1.0 + (QEZ4 + QEZ5 * GAMZ) * (2.0 / np.pi) * np.arctan(BT * CT * ALPHAT)))

    BR = QBZ9 * L['LKY'] / L['LMUY'] + QBZ10 * BY * CY
    DR = (FZ * ((QDZ6 + QDZ7 * DFZ) * L['LRES'] + (QDZ8 + QDZ9 * DFZ) * GAMZ) * r0 * L['LMUY'])

    TRAIL = DT * np.cos(CT * np.arctan(BT * ALPHAT - ET * (BT * ALPHAT - np.arctan(BT * ALPHAT)))) * np.cos(ALPHA)
    MZR = DR * np.cos(np.arctan(BR * ALPHAR)) * np.cos(ALPHA)

    return -TRAIL * FY0 + MZR

def FX(X, p, fz0):
    SR, FZ, IA = X
    KAPPA = SR
    GAMMA = np.radians(IA)

    (PCX1, PDX1, PDX2, PDX3, PEX1, PEX2, PEX3, PEX4,
     PKX1, PKX2, PKX3, PHX1, PHX2, PVX1, PVX2) = p

    FZ0 = fz0 * L['LFZO']
    DFZ = (FZ - FZ0) / FZ0
    GAMX = GAMMA * L['LGAX']

    SHX = (PHX1 + PHX2 * DFZ) * L['LHX']
    KAPPAX = KAPPA + SHX
    SVX = FZ * (PVX1 + PVX2 * DFZ) * L['LVX']
    CX = PCX1 * L['LCX']
    MUX = (PDX1 + PDX2 * DFZ) * (1.0 - PDX3 * GAMX ** 2) * L['LMUX']
    DX = MUX * FZ
    KX = FZ * (PKX1 + PKX2 * DFZ) * np.exp(PKX3 * DFZ) * L['LKX']
    BX = KX / (CX * DX + 1e-9)
    EX = (PEX1 + PEX2 * DFZ + PEX3 * DFZ ** 2) * (1.0 - PEX4 * np.sign(KAPPAX)) * L['LEX']

    return DX * np.sin(CX * np.arctan(BX * KAPPAX - EX * (BX * KAPPAX - np.arctan(BX * KAPPAX)))) + SVX
