def getPdfWeight(fname):
  with open('pdfuncerts_v3.log') as f:
    for l in f.readlines():
      content0, content1, content2 = l.split(' ')[0], l.split(' ')[1], l.split(' ')[2]
      if fname in content0:
        return content1, content2

