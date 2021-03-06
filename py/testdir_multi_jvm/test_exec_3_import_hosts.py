import unittest
import random, sys, time, os
sys.path.extend(['.','..','py'])
import h2o, h2o_cmd, h2o_hosts, h2o_browse as h2b, h2o_import as h2i, h2o_exec as h2e

zeroList = []
for i in range(8):
    zeroList.append('ColumnRes' + str(i) + ' = 0')
    zeroList.append('MatrixRes' + str(i) + ' = 0')
    zeroList.append('ScalarRes' + str(i) + ' = 0')

# FIX! put these in 3?
# 'randomBitVector' ?? hardwire size to 19?
# 'randomFilter'
# 'log"
# do we have to restrict ourselves?
# 'factor' (hardware the enum colum to col 53
# bug?
#        'Result<n> = randomFilter(<keyX>[<col1>],<row>)',
#        'MatrixRes<n> = slice(<keyX>[<col1>],<row>)',
exprList = [
        'ColumnRes<n> = factor(<keyX>[53]) + ColumnRes<n-1>',
        'ColumnRes<n> = randomBitVector(19,0,1) + ColumnRes<n-1>',
        'ColumnRes<n> = log(<keyX>[<col1>]) + ColumnRes<n-1>',
        'ColumnRes<n> = <keyX>[<col1>] + <keyX>[<col2>] + <keyX>[2]',
        'MatrixRes<n> = colSwap(<keyX>,<col1>,(<keyX>[2]==0 ? 54321 : 54321))',
        'ColumnRes<n> = <keyX>[<col1>]',
        'ScalarRes<n> = min(<keyX>[<col1>])',
        'ScalarRes<n> = max(<keyX>[<col1>]) + ColumnRes<n-1>',
        'ScalarRes<n> = mean(<keyX>[<col1>]) + ColumnRes<n-1>',
        'ScalarRes<n> = sum(<keyX>[<col1>]) + ScalarRes<n-1>',
    ]

class Basic(unittest.TestCase):
    def tearDown(self):
        h2o.check_sandbox_for_errors()

    @classmethod
    def setUpClass(cls):
        global SEED, localhost
        SEED = h2o.setup_random_seed()
        localhost = h2o.decide_if_localhost()
        if (localhost):
            h2o.build_cloud(3,java_heap_GB=4)
        else:
            h2o_hosts.build_cloud_with_hosts()

    @classmethod
    def tearDownClass(cls):
        # wait while I inspect things
        # time.sleep(1500)
        h2o.tear_down_cloud()

    def test_exec_import_hosts(self):
        # make the timeout variable per dataset. it can be 10 secs for covtype 20x (col key creation)
        # so probably 10x that for covtype200
        if localhost:
            maxTrials = 200
            csvFilenameAll = [
                ("covtype.data", "cA", 15),
            ]
        else:
            maxTrials = 20
            csvFilenameAll = [
                ("covtype.data", "cB", 15),
                ("covtype20x.data", "cD", 60),
            ]

        ### csvFilenameList = random.sample(csvFilenameAll,1)
        csvFilenameList = csvFilenameAll
        h2b.browseTheCloud()
        lenNodes = len(h2o.nodes)
        importFolderPath = "standard"
        cnum = 0
        for (csvFilename, hex_key, timeoutSecs) in csvFilenameList:
            cnum += 1
            # creates csvFilename.hex from file in importFolder dir 
            csvPathname = importFolderPath + "/" + csvFilename
            parseResult = h2i.import_parse(bucket='home-0xdiag-datasets', path=csvPathname, 
                hex_key=hex_key, timeoutSecs=2000)
            print csvFilename, 'parse time:', parseResult['response']['time']
            print "Parse result['destination_key']:", parseResult['destination_key']

            # We should be able to see the parse result?
            inspect = h2o_cmd.runInspect(None, parseResult['destination_key'])

            print "\n" + csvFilename
            h2e.exec_zero_list(zeroList)
            # we use colX+1 so keep it to 53
            # we use factor in this test...so timeout has to be bigger!
            h2e.exec_expr_list_rand(lenNodes, exprList, hex_key, 
                maxCol=53, maxRow=400000, maxTrials=maxTrials, timeoutSecs=(timeoutSecs))


if __name__ == '__main__':
    h2o.unit_main()
