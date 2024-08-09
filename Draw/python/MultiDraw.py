# MultiDraw.py

import ROOT
import re
import string
import random
from array import array
import numpy as np
import os

file_directory = os.path.dirname(__file__)
lib_path = os.path.join(os.path.dirname(file_directory), 'lib/libMultiDraw.so')
ROOT.gSystem.Load(lib_path)
ROOT.gInterpreter.Declare('''
        extern void MultiDraw(TTree *inTree,
            TObjArray *Formulae, TObjArray *Weights, TObjArray *Hists,
            UInt_t ListLen);
''')


def split_vals(vals):
    """Converts a string '1:3|1,4,5' into a list [1, 2, 3, 4, 5]"""
    res = set()
    first = vals.split(',')
    for f in first:
        second = re.split('[:|]', f)
        # print second
        if len(second) == 1:
            res.add(float(second[0]))
        if len(second) == 3:
            x1 = float(second[0])
            while x1 < float(second[1]) + 1E-5:
                res.add(x1)
                x1 += float(second[2])
    return sorted([x for x in res])


def randomword(length):
    return ''.join(random.choice(string.lowercase) for i in range(length))


def MakeTObjArray(theList, takeOwnership=True):
    """Turn a python iterable into a ROOT TObjArray"""

    result = ROOT.TObjArray()
    if takeOwnership:
        result.SetOwner()

    # Make PyROOT give up ownership of the things that are being placed in the
    # TObjArary. They get deleted because of result.SetOwner()
    for item in theList:
        if takeOwnership:
            ROOT.SetOwnership(item, False)
        result.Add(item)

    return result


def GetBinningArgs(arg, is_variable):
    if is_variable:
        binning = split_vals(arg)
        return [len(binning) - 1, array('d', binning)]
    else:
        str_binning = [x.strip() for x in arg.split(',')]
        binning = []
        if len(str_binning) == 3:
            step = (float(str_binning[2]) - float(str_binning[1])) / float(str_binning[0])
            binning = np.arange(float(str_binning[1]), float(str_binning[2]) + step / 2,step)
            return [len(binning) - 1, array('d',binning)]
        else:
            return None


def MultiDraw(self, Formulae, Compiled=False):
    results, formulae, weights, formulaeStr, weightsStr = [], [], [], [], []

    # lastFormula, lastWeight = None, None

    for i, origFormula in enumerate(Formulae):

        # Expand out origFormula and weight, otherwise just use weight of 1.
        if type(origFormula) is tuple:
            origFormula, weight = origFormula
        else:
            origFormula, weight = origFormula, "1"

        # Our way is simpler, require each variable to end in (...) or [...] to give the binning
        # and always create a new hist

        split_var = origFormula.split(';')
        origFormula = split_var[0]
        print("Formula: ", origFormula, weight)

        var_binned_x = False
        var_binned_y = False
        var_binned_z = False
        is_2d = False
        is_3d = False

        if origFormula[-1] == ')':
            pos_open = origFormula.rfind('(')
            pos_close = origFormula.rfind(')')
        if origFormula[-1] == ']':
            var_binned_x = True
            pos_open = origFormula.rfind('[')
            pos_close = origFormula.rfind(']')
        if pos_open == -1 or pos_close == -1 or pos_open > pos_close:
            raise RuntimeError('You bus')
        bin_args_x = GetBinningArgs(origFormula[pos_open + 1:pos_close], var_binned_x)
        formula = origFormula[:pos_open].strip()

        # Check if this is a 2D histogram with syntax
        # [var_y],[var_x][binning_y],[binning_x]
        if formula[-1] == ',':
            is_2d = True
            if formula[-2] == ')':
                pos_open_y = formula.rfind('(')
                pos_close_y = formula.rfind(')')
            if formula[-2] == ']':
                var_binned_y = True
                pos_open_y = formula.rfind('[')
                pos_close_y = formula.rfind(']')
            if pos_open_y == -1 or pos_close_y == -1 or pos_open_y > pos_close_y:
                raise RuntimeError('You bus')
            bin_args_y = GetBinningArgs(formula[pos_open_y + 1:pos_close_y], var_binned_y)
            formula = origFormula[:pos_open_y].strip()
            # Check if this is a 3D histogram with syntax
            # [var_y],[var_x],[var_z][binning_y],[binning_x],[binning_z]
            if formula[-1] == ',':
                is_3d = True
                if formula[-2] == ')':
                    pos_open_z = formula.rfind('(')
                    pos_close_z = formula.rfind(')')
                if formula[-2] == ']':
                    var_binned_z = True
                    pos_open_z = formula.rfind('[')
                    pos_close_z = formula.rfind(']')
                if pos_open_z == -1 or pos_close_z == -1 or pos_open_z > pos_close_z:
                    raise RuntimeError('You bus')
                bin_args_z = GetBinningArgs(formula[pos_open_z + 1:pos_close_z], var_binned_z)
                formula = formula[:pos_open_z].split(',')
            else:
                formula = formula[:pos_open_y].split(',')
        else:
            formula = [formula]

        ROOT.TH1.AddDirectory(False)
        if not is_2d and not is_3d:
            hist = ROOT.TH1D(origFormula + ':' + weight, origFormula, *bin_args_x)
        elif not is_3d:
            hist = ROOT.TH2F(origFormula + ':' + weight, origFormula, *(bin_args_x + bin_args_y))
        else:
            hist = ROOT.TH3F(origFormula + ':' + weight, origFormula, *(bin_args_x + bin_args_y + bin_args_z))

        if len(split_var) > 1:
            hist.GetXaxis().SetTitle(split_var[1])
        if len(split_var) > 2:
            hist.GetXaxis().SetTitle(split_var[2])
            hist.GetYaxis().SetTitle(split_var[1])
        if len(split_var) > 2:
            hist.GetXaxis().SetTitle(split_var[3])
            hist.GetYaxis().SetTitle(split_var[2])
            hist.GetZaxis().SetTitle(split_var[1])

        if is_2d:
            results.append(ROOT.TObject())
        if is_3d:
            results.append(ROOT.TObject())
        results.append(hist)

        # The following two 'if' clauses check that the next formula is different
        # to the previous one. If it is not, we add an ordinary TObject.
        # Then, the dynamic cast in MultiDraw.cxx fails, giving 'NULL', and
        # The previous value is used. This saves the recomputing of identical
        # values

        for form in formula:
            f = ROOT.TTreeFormula("formula%i" % i, form, self)
            f.SetTitle(form)
            if not f.GetTree():
                raise RuntimeError("TTreeFormula didn't compile: " + form)
            f.SetQuickLoad(True)
            formulae.append(f)
            formulaeStr.append(form)

            f = ROOT.TTreeFormula("weight%i" % i, weight, self)
            f.SetTitle(weight)
            if not f.GetTree():
                raise RuntimeError("TTreeFormula didn't compile: " + weight)
            f.SetQuickLoad(True)
            weights.append(f)
            weightsStr.append(weight)

    from ROOT import MultiDraw as _MultiDraw
    from time import time
    start = time()

    # Ensure that formulae are told when tree changes
    fManager = ROOT.TTreeFormulaManager()
    for formula in formulae + weights:
        if type(formula) is ROOT.TTreeFormula:
            fManager.Add(formula)

    fManager.Sync()
    self.SetNotify(fManager)

    # Draw everything!
    _MultiDraw(self,
               MakeTObjArray(formulae),
               MakeTObjArray(weights),
               MakeTObjArray(results, takeOwnership=False),
               len(formulae))

    print("Took %.2fs" % (time() - start), " " * 20)
    return results


ROOT.TTree.MultiDraw = MultiDraw
