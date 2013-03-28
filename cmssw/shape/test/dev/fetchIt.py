#!/usr/bin/python

import re, urllib, urllib2
import xlrd
import csv
import tempfile
import os

class Spreadsheet(object):
    def __init__(self, key):
        super(Spreadsheet, self).__init__()
        self.key = key

class Client(object):
    def __init__(self, email, password):
        super(Client, self).__init__()
        self.email = email
        self.password = password

    def _get_auth_token(self, email, password, source, service):
        url = "https://www.google.com/accounts/ClientLogin"
        params = {
            "Email": email, "Passwd": password,
            "service": service,
            "accountType": "HOSTED_OR_GOOGLE",
            "source": source
        }
        req = urllib2.Request(url, urllib.urlencode(params))
        return re.findall(r"Auth=(.*)", urllib2.urlopen(req).read())[0]

    def get_auth_token(self):
        source = type(self).__name__
        return self._get_auth_token(self.email, self.password, source, service="wise")

    def download(self, spreadsheet, gid=0, format="csv"):
        url_format = "https://spreadsheets.google.com/feeds/download/spreadsheets/Export?key=%s&exportFormat=%s&gid=%i"
        headers = {
            "Authorization": "GoogleLogin auth=" + self.get_auth_token(),
            "GData-Version": "3.0"
        }
        req = urllib2.Request(url_format % (spreadsheet.key, format, gid), headers=headers)
        return urllib2.urlopen(req)

    def download_xls(self, spreadsheet, format="xls"):
        url_format = "https://spreadsheets.google.com/feeds/download/spreadsheets/Export?key=%s&exportFormat=%s"
        headers = {
            "Authorization": "GoogleLogin auth=" + self.get_auth_token(),
            "GData-Version": "3.0"
        }
        req = urllib2.Request(url_format % (spreadsheet.key, format), headers=headers)
        return urllib2.urlopen(req)

if __name__ == "__main__":
    import getpass
    import csv

    tempname = 'dummy.xls'

#     tmp_xls = tempfile.NamedTemporaryFile()

    ssname = 'Test1'
    email = "alessandro.thea@gmail.com" # (your email here)
    password = getpass.getpass()
    spreadsheet_id = "0AmbqMj_rTADpdDNhVnZfckxXLTZYcGNfVGZaMTRSQ3c" # (spreadsheet id here)

    # Create client and spreadsheet objects
    gs = Client(email, password)
    ss = Spreadsheet(spreadsheet_id)

    print "downloading...",
    # Request a file-like object containing the spreadsheet's contents
    xls_file = gs.download_xls(ss)

#     tmp_xls.write(xls_file.read())

    f = open(tempname, 'wb')
    f.write(xls_file.read())
    f.close()
    print "done"
#     print "converting",tmp_xls.name,'to',ssname+'csv...',
    print "converting",tempname,'to',ssname+'csv...',

#     wb = xlrd.open_workbook(tmp_xls.name)
    wb = xlrd.open_workbook(tempname)
    sh = wb.sheet_by_name(ssname)
    your_csv_file = open(ssname+'.csv', 'wb')
    wr = csv.writer(your_csv_file, quoting=csv.QUOTE_ALL)

    for rownum in xrange(sh.nrows):
        wr.writerow(sh.row_values(rownum))

    your_csv_file.close()
    print "done"
#     tmp_xls.close()

    os.unlink(tempname)



