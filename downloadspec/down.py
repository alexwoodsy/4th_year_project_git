import wget

f = open("downloadtest_url.txt", "r")
fitsurl = f.read().split()
fitsurl = tuple(fitsurl)
f.close()

for i in range(0,len(fitsurl)):
    urldown = fitsurl[i]
    wget.download(urldown, "downloads")

#testing commit feature
