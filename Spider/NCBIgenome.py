import bs4
import csv
import requests

i = input("请输入关键词：")
url = "https://www.ncbi.nlm.nih.gov/genome?term="+str(i)+"&EntrezSystem2.PEntrez.Genome2.Genome2_ResultsPanel.Genome2_DisplayBar.PageSize=200"
html = requests.get(url, timeout = 100).text

soup = bs4.BeautifulSoup(html, "html.parser")
urllist = soup.select('.title a')
for u in urllist:
    print(u.get('href'))
    content = requests.get('https://www.ncbi.nlm.nih.gov' + u.get('href'), timeout = 100).text
    soup1 = bs4.BeautifulSoup(content, "html.parser")
    GenomeTitle = soup1.select('.GenomeTitle')[0].get_text()
    ftpurl = soup1.select('.shifted a')[0].get('href')
    Assembly = ''
    try:
        for x in range(3):
            if soup1.select('.summary table tr')[x].get_text().find('Assembly:') > -1:
                Assembly = soup1.select('.summary table tr')[x].get_text()
                break
    except Exception as e:
        raise e
    print(Assembly)
    Statistics = ''
    Statistics1 = ''
    Statistics2 = ''
    try:
        for x in range(7):
            if soup1.select('.summary table tr')[x].get_text().find('Statistics') > -1:
                Statistics = soup1.select('.summary table tr')[x].get_text()
                if soup1.select('.summary table tr')[x+1].get_text().find('count') > -1:
                    Statistics1 = soup1.select('.summary table tr')[x+1].get_text()
                try:
                    if soup1.select('.summary table tr')[x+2].get_text().find('GC') > -1:
                        Statistics2 = soup1.select('.summary table tr')[x+2].get_text()
                except Exception as e:
                    raise e
                break
    except Exception as e:
        raise e
    print(Statistics+Statistics1+Statistics2)
    result = []
    result.append(GenomeTitle)
    result.append(Assembly.replace(u'\xa0', u' '))
    result.append((Statistics+Statistics1+Statistics2).replace(u'\xa0', u' '))
    result.append(ftpurl)
    result.append(" ")
    with open("result.csv", "a", newline='', encoding = 'utf-8') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(result)
print(GenomeTitle,Assembly,Statistics,ftpurl)

