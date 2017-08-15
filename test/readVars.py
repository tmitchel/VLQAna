def getPdfWeight(fname, channel):
  with open('pdfuncerts_'+channel+'.log') as f:
    for l in f.readlines():
      content0, content1, content2, content3, content4 = l.split(' ')[0], l.split(' ')[1], l.split(' ')[2], l.split(' ')[3], l.split(' ')[4]
      if fname in content0:
        return content1, content2, content3, content4
    else:
      return 1, 1, 1, 1

