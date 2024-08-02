#!/usr/bin/python
import sys 
import os
import numpy as np
import re

dicto={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X'}

motifs_list={'TAV',
'QQQ',
'VLH',
'TTP',
'AFG',
'KPA',
'WGQ',
'GLY',
'LFG',
'GGY',
'GGG',
'EDD',
'VPP',
'PRG',
'NQY',
'PGT',
'FGS',
'QPP',
'FGA',
'TPG',
'TPK',
'TPS',
'TPQ',
'TPP',
'KAK',
'APP',
'QYN',
'SFG',
'RTP',
'DLY',
'AFS',
'SSY',
'SSS',
'SSH',
'HGG',
'QTP',
'DGI',
'PDY',
'QFQ',
'HVP',
'SKC',
'AAP',
'GYS',
'GYQ',
'RPK',
'FKG',
'NRG',
'QAT',
'TQT',
'PAA',
'PAF',
'SYG',
'PHG',
'RGR',
'SKP',
'RGA',
'EMP',
'QPH',
'NVP',
'PSV',
'DVY',
'DVP',
'KDQ',
'NDE',
'SRR',
'QNY',
'TDD',
'PPA',
'QRR',
'GSY',
'DRK',
'LYQ',
'YDV',
'KVH',
'EEE',
'EED',
'AMA',
'GGW',
'KCG',
'YQG',
'YQQ',
'RRL',
'TEA',
'CGS',
'HHS',
'DNI',
'LEP',
'YNN',
'YNQ',
'DAL',
'YGG'
}
text_file1= open("motifs1.txt", "w")

clv={'[DSTE][^P][^DEWHFYC]D[GSAN]',
'[ILV]..R[VF][GS].',
'(.RK)|(RR[^KR])',
'R.[RK]R.',
'[KR]R.',
'KR.',
'R...[KR]R.',
'[RK].[AILMFV][LTKF].',
'S[IVLMH]E[IVPFMLYAQR]GR.',
'E[IMPVL][MLVP]R.',
'Q[MLVI]DG..[DE]',
}
text_file2= open("motifsclv.txt", "w")
doc={'F..[FWY][ST][FY]',
'F..[FWY][DE][FY]',
'F..F$',
'.R..[PGAV][DEIP]G.',
'[MPVLIFWYQ].(T)P..',
'[RK].L.{0,1}[FYLIVMP]',
'V[ED]P[^P][RK]FA[^P]ELI[^P]RLE[^P][VIL]',
'[RK].{2,4}[LIVP]P.[LIV].[LIVMF]|[RK].{2,4}[LIVP].P[LIV].[LIVMF]',
'F.[FY]P',
'[KR]{0,2}[KR].{0,2}[KR].{2,4}[ILVM].[ILVF]',
'([LIV][^P][^P][RK]....[LIVMP].[LIV].[LIVMF])|([LIV][^P][^P][RK][RK]G.{4,7}[LIVMP].[LIV].[LIVMF])',
'[RK]P[^P][^P]L.[LIVMF]',
'[RK].{2,4}[LIVMP].[LIV].[LIVMF]',
'[RK][^P][^P][LIM].L.[LIVMF].',
'[LIVMPFA].[LIV].{1,2}[LIVMP].{4,6}[LIV]..[RK][RK]',
'[DEN][DEN].{2,3}[ILMVA][DEN][DEN]L',
'R[^P][DEQ]Q[VIL]([RK][^P]|[^P][RK])[YW]',
'..[RK].{0,1}[VIL][^P][FW].',
'.[GS]IL[KR][^DE]',
'([LMFYWIC]..I.E)|(L..[IVLWC].E).',
'L.[LIVAPM]P',
'.P[^P]I[^P][IV][^P]',
'RF[^P][IV].',
'[PA][^P][^FYWIL]S[^P]',
'P.E[^P].S[^P]',
'K...K',
'F[EDQS][MILV][ED][MILV]((.{0,1}[ED])|($))',
'...([ST])P.'
}
text_file3= open("motifsdoc.txt", "w")
deg={
'.R..L..[LIVM].',
'.KEN.',
'.[ILM]R$',
'[STDE]{1,3}.{0,2}[TSDE].{2,3}VP[STDE]G{0,1}[FLIMVYPA]',
'[NQ]{0,1}..[ILMV][ST][DEN][FY][FY].{2,3}[KR]{2,3}[^DE]',
'[NQ]{0,1}..[ILMV]T[DEN][HMFY][FMY].{2,3}[KR]{2,3}[^DE]',
'[AP]P[MV][IM]V',
'[DNS].[DES][TNS]GE',
'QD.DLGV',
'E.EE.E[AV]DQH',
'F[^P]{3}W[^P]{2,3}[VIL]',
'^M{0,1}[FYLIW][^P]',
'^M{0,1}[RK][^P].',
'^M{0,1}([ED]).',
'^M{0,1}([NQ]).',
'^M{0,1}(C).',
'[IL]A(P).{6,8}[FLIVM].[FLIVM]',
'..[RK][RK].SL..F[FLM].[RK]R[HRK].[RK].',
'[LIVMP].{0,2}(T)P..([ST])',
'[LIVMP].{0,2}(T)P..E',
'..[DE].(T)P.K',
'.[VLIA][VLI]GWPP[VLI]...R.',
'D(S)G.{2,3}([ST])',
'.P.A.V.P[^P]',
'[AVP].[ST][ST][ST]'
}

text_file4= open("motifsdeg.txt", "w")


lig={'R[^DE]{0,2}[^DEPG]([ST])(([FWYLMV].)|([^PRIKGN]P)|([^PRIKGN].{2,4}[VILMFWYP]))',
'R[^DE]{0,2}[^DEPG]([ST])[^P]{0,1}$',
'[IL]..[^P][^P][^P][^P]R.....[IL]..[^P][^P][ILV][ILM]',
'R..[ILVMF][ILMVF][^P][^P][ILVM].{4,7}L(([KR].)|(NK))[VATI]',
'[^R]..((.[ILMVF])|([ILMVF].))[^P][^P][ILVM].{4,7}L(([KR].)|(NK))[VATIGS]',
'P.LP.[IL].{1,3}[VLF]',
'F.D.F',
'DP[FW]',
'[ILVMF].[ILMVP][FHY].[DE]',
'[KR]..[ILVM][FHY].[DE]',
'[DE]R[YFH][ILFVM][PAG].R',
'DR[YFH][ILFVM][PA]..',
'[DE][DES][DEGAS]F[SGAD][DEAP][LVIMFD]',
'....[LIFVYMTE][ASGC][^P]{2}L[^P]{2}[IVMTL][GACS][D][^P][FVLMI].',
'^M{0,1}[AS]...',
'^M{0,1}A.P.',
'DA.P.',
'^M{0,1}A.[AP].',
'DA.G.',
'.(S)..F',
'.(S)..F.K',
'.(S)..Y$',
'[ACLIVTM][^P][^P][ILVMFCT]Q[^P][^P][^P][RK][^P]{4,5}[RKQ][^P][^P]',
'((SP)|([ED].{0,1}))[IV]W[IVL].R',
'[ED].{0,2}[ED].{0,2}[EDQ].{0,1}[YF]$',
'.W[RK][DE]GCY$',
'[DE][DEN][DEN]D[GDN]Y.P..',
'L[IVLMF].[IVLMF][DE]',
'.[NP]W[DES].W',
'[FY][^P].[WFY][^P]DY..L',
'L[^P]{2,2}[HI]I[^P]{2,2}[IAV][IL]',
'EP[IL]Y[TAG]',
'[AFILMPTVW]W[FHILMPSTVW]P',
'(P[LVIPME][DENS][LM][VASTRG])|(G[LVIPME][DENS][LM][VASTRG]((K)|(.[KR])))',
'^M[MIL].[MIL]',
'[^P].[KR].TQT',
'.A.GPP.{2,3}Y.',
'P[PG]{0,1}YP.{1,6}Y[QS]{0,1}P',
'P.P.{0,1}GF',
'.NPF.',
'.[FYH].[IVM][^WFYP][^WFYP][ILM][ILMV].',
'Y....L[VILMF]',
'Y.PP.[ILMV]R',
'([FYWL]P.PP)|([FYWL]PP[ALIVTFY]P)',
'PP..F',
'[FY].[FW].....[LMVIF]P.P[DE]',
'[LV][DE][^P][LM][LM][^P][^P]L[^P]',
'..(T)..[ILV].',
'..(T)..[DE].',
'W.[VIL].[ST].KA{0,1}T...W',
'[FYLIMV].FG[DES]F',
'[ILMV]...[ILMVF]..[ILMVA][ILMVA].[KR]R..[ILMVA]',
'[EN][FYLW][NSQ].EE[ILMVF][^P][LIVMFA]',
'(([CP]PP)|(PP[TP]))[ST]P[^P][TS]{0,1}',
'[QHR].{0,1}P[PL]PP[GS]H[RH]',
'[DE]H.Y',
'[FY][DEP]WM',
'P[MVLIRWY]V[MVLIAS][LM]',
'G[FL]PGER..G',
'NGR',
'..L.I(S)',
'[VILMFT]K.EP.[DE]',
'[VILMFT]K.EP.{2,3}[DE]',
'[VILMFT]K.EP....[DE]',
'[LMTAFSRI][^KRG]W[DE].{3,5}[LIVMFPA]',
'[EDST].{0,2}[WFY]..P',
'[EDST].{0,2}[WFY]..[ILV]',
'[EDST].{0,2}LVV',
'[EDST].{0,2}[WFY]..[ILVFY]',
'([VILA]..N.I[RK])|([VILA].PN.IG.{0,6}[RK])',
'[LM]YP...[LI][^P][^P][LI]',
'[LM]YP.[LI]',
'[KR][IV][LV].....P',
'GRYFG',
'[LFVAIMW].{3,5}[DE][FY][IL][SAPGK][FL].{3,6}[DE]{3}',
'Q[RGQ]DF[LI][PS]L[DE]',
'P.L.P',
'PP.LI',
'[LMV]P.LE',
'F..A[ILV]..A..[ILV]',
'[^P]L[^P][^P]LL[^P]',
'.F[^P][^P][KRIL]H[^P][^P][YLMFH][^P]...',
'....WF..L',
'..[LFP][NS][PIVTAFL].A..(([FY].[PYLF])|(W..)).',
'((WPP)|([FL][PV][APQ]))EF.PG.PWKG.',
'((^.{0,3})|(Q)).[^FHWY][ILM][^P][^FHILVWYP][HFM][FMY]..',
'...[ST].[ACVILF]$',
'...[VLIFY].[ACVILF]$',
'...[DE].[ACVILF]$',
'W...[FY]',
'F...[WF]',
'LV.EF[LM]',
'MM[NDE][EDNAG]F[LMA]',
'L..LL...L..F',
'.P[TS]AP.',
'(.[^P].NP.[FY].)|(.[ILVMFY].N..[FY].)',
'(.[^P].NP.(Y))|(.[ILVMFY].N..(Y))',
'([DEST]|^).{0,4}[LI].C.E.{1,4}[FLMIVAWPHY].{0,8}([DEST]|$)',
'..[LIMV]..[LM][FY]D.',
'RGD',
'[^P][MIALVF][^P][^P][NSHRA]R[^P][^P][ASV][^P][^P][RKLIA][RQLIVA]',
'[^P][MIAL][^P][^P][NKS][KRLQH][^P][^P]A[^P][^P][RKLI][RKL][^P][^P][KR]',
'R[MIVAS][^P][^P][NQ][KRL][^P][^P]A[^P][^P][RK]',
'[KRS]I[^P][^P][NK][KR][^P][^P]A[^P][^P][RKL][RKL][^P][^P][RK]',
'.[ILVM]LG..P.',
'(Y).N.',
'(Y)[IV].[VILP]',
'(Y)[QDEVAIL][DENPYHI][IPVGAHS]',
'(Y)..Q',
'(Y)[VLTFIC]..',
'G(Y)[KQ].F',
'[RKY]..P..P',
'P..P.[KR]',
'...[PV]..P',
'KP..[QK]...',
'P..DY',
'[LIV]..[LM]L.AA.[FY][LI]',
'[FHYM].A[AV].[VAC]L[MV].[MI]',
'[FA].[LA][LV][LVI]..[AM]',
'[ED][LIV]NNN[^P]',
'[SV][CY]GH[LIF][LAST][GAIV].',
'[DEST]{1,10}.{0,1}[VIL][DESTVILMA][VIL][VILM].[DEST]{0,5}',
'[DEST]{0,5}.[VILPTM][VIL][DESTVILMA][VIL].{0,1}[DEST]{1,10}',
'([KR][^ED]{0,5}[ST].IP[^ED]{5,5})|([^ED]{5,5}[ST].IP[^ED]{0,5}[KR])',
'EEVD$',
'[PSAT].[QE]E',
'P.Q..D',
'..P.E..[FYWHDE].',
'[FY].L.P',
'[DEN]..(Y)..[LI].{6,12}(Y)..[LI]',
'[ILV].(Y)..[ILV]',
'..T.(Y)..[IV]',
'[ILM][ILMF].{1,2}[ILM].{0,4}K',
'[ND].WGI.[LIV][VMLI].{0,1}[ED]',
'[KR]{1,4}[KR].[KR]W.',
'[VMIL][MILFYPA][^P][TASKHC][AVSC][^P][^P][ILVM][^P][^P][^P][LMTVI][^P][^P][LMV][ILVMA][^P][^P][AIVLMT]',
'[ED].{0,3}[VIL]D[VI]',
'[EDSTY].{0,4}[VIPLA][TSDEKR][ILVA]',
'[SCA]AR[STCA][EQR][PGILVM][HYFQNKRLVI]',
'[STCA][CSAGV]R[STCAV][EQR][PGALV][LFYHRK]',
'[SCA][AFWHSV][KR][TAS][DEQR][GP][RKYFWIVAM]..[IVM]',
'ES[RK][FY].F[HR][PST][IVLM][DES][DE]',
'[WFY]RP[WFY].{0,7}$',
'[WFY][KR]P[WFY]',
'PP.Y',
'PPLP',
}
text_file5= open("motifslig.txt", "w")
trg={'[DE].{1,2}F[^P][^P][FL][^P][^P][^P]R',
'QV.P.$',
'RV.P.',
'Y..[LMVIF]',
'([LIVMFYWPR]R[^YFWDE]{0,1}R)|(R[^YFWDE]{0,1}R[LIVMFYWPR])',
'K.{0,1}K.{2,3}$',
'[DE].{0,4}E[FY][FYK]D[AC].[ESTD]',
'[KRHQSAP][DENQT]EL$',
'Q.{6,6}FF.{6,7}$',
'[DERQ]...L[LVI]',
'[DET]E[RK].PL[LI]',
'D..LL.{1,2}$',
'S[LW]LD[DE]EL[LM]',
'([DEQ].{0,1}[LIM].{2,3}[LIVMF][^P]{2,3}[LMVF].[LMIV].{0,3}[DE])|([DE].{0,1}[LIM].{2,3}[LIVMF][^P]{2,3}[LMVF].[LMIV].{0,3}[DEQ])',
'[KR][KR].{7,15}[^DE]((K[RK])|(RK))(([^DE][KR])|([KR][^DE]))[^DE]',
'[^DE]((K[RK])|(RK))[KRP][KR][^DE]',
'[^DE]((K[RK])|(RK))(([^DE][KR])|([KR][^DE]))(([PKR])|([^DE][DE]))',
'(([PKR].{0,1}[^DE])|([PKR]))((K[RK])|(RK))(([^DE][KR])|([KR][^DE]))[^DE]',
'(.[SAPTC][KRH][LMFI]$)|([KRH][SAPTC][NTS][LMFI]$)',
'^.{1,40}R[^P][^P][^P][LIV][^P][^P][HQ][LIF]',
}
text_file6= open("motifstrg.txt", "w")
mod={'C.([DN]).{4,4}[FY].C.C',
'(C)[^DENQ][LIVM].$',
'...([ST])P[RK]',
'...([ST])P.[KR]',
'...([ST])P..[RK]',
'S..([ST])...',
'...([ST])..E',
'(W)..W',
'(.)G[RK][RK]',
'[ED]{0,3}.(S)[GA].',
'...([ST])...[ST]',
'H.[KR]..([ST])[^P]',
'[FLM][^P][^P]([ST])[^DEP][^DE]',
'[WYPCAG][^P][^P]([ST])[IFCVML][KRHYF]',
'.(N)[^P][ST]..',
'(N)[^P]C',
'^M{0,1}(G)[^EDRKHPFYW]..[STAGCN][^P]',
'C.{3,5}([ST])C',
'C.(S).PC',
'...([ST])Q..',
'[RK]..(S)[VI]..',
'[RK][RK].([ST])[^P]..',
'.R.([ST])[^P]..',
'R.R..([ST])[^P]..',
'.[DE].([ST])[ILFWMVA]..',
'...([ST])P..',
'G(C)M[GS][CL][KP]C',
'^M{0,1}G(C)..S[AKS]',
'[VILMAFP](K).E',
'[SDE].{0,5}[DE].(K).{0,1}[AIFLMPSTV]',
'[TAD][EA].Q(Y)[QE].[GQA][PEDLS]',
'..[RKTC][IVL]Y[TQHS](Y)[IL]QSR',
'[ETA](C)[QERK]..F...RWNC[ST]'
}
text_file7= open("motifsmod.txt", "w")


def motifslist(seq,list2):
      motifs=[]
      G=[]
      G=[0] * (len(seq))
      for i in range(len(seq)-4):
         p= str(seq[i])+str(seq[i+1])+str(seq[i+2])
         if p in list2:
              G[i]+=1
              G[i+1]+=1
              G[i+2]+=1
      motifs=G
      return motifs

def findmotifs(seq,list2):
        motifs=[]
        G=[]
        G=[0] * (len(seq))
        for t in list2:
         regex=re.compile(t)
         it = re.finditer(regex, seq)
         for match in it:
              s=match.group()
              x=seq.find(s)
              if x>0:
               for i in range(x,x+len(match.group())-1):
                 G[i]+=1
        return G 
                  
          

name2=str(sys.argv[1])
fp=open(name2,'r')


for line in fp:
      line=line.replace('\n','')
      line=line.replace(' ','')
      finnok=motifslist(line,motifs_list)
      for i in finnok:
           text_file1.write("%s \n" %(i))
      text_file1.close()
      finnok=findmotifs(line,clv)
      for i in finnok:
           text_file2.write("%s \n" %(i))
      text_file2.close()
      finnok=findmotifs(line,doc)
      for i in finnok:
           text_file3.write("%s \n" %(i))
      text_file3.close()
      finnok=findmotifs(line,deg)
      for i in finnok:
           text_file4.write("%s \n" %(i))
      text_file4.close()
      finnok=findmotifs(line,lig)
      for i in finnok:
           text_file5.write("%s \n" %(i))
      text_file5.close()
      finnok=findmotifs(line,trg)
      for i in finnok:
           text_file6.write("%s \n" %(i))
      text_file6.close()
      finnok=findmotifs(line,mod)
      for i in finnok:
           text_file7.write("%s \n" %(i))
      text_file7.close()
      
