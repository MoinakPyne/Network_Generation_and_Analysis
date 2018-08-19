import urllib
import unicodedata
import zen
import time
import itertools
import re
from bs4 import BeautifulSoup

# page content starts with: 
PAGE_START = '<ul class="results">'
PAGE_END = '</ul>'
PAPER_START = '<li aria-describedby='
PAPER_END = '</li>'
CONTENT_START = '<div class="txt">'
CONTENT_END = '</div>'
TITLE_START = '<h3>'
TITLE_END = '</h3>'
TITLE_CONTENT_START = '<a'
TITLE_CONTENT_END = '</a>'
PAPER_SPAN_START = '<span id="art-abs-title'
PAPER_SPAN_END = '</span>'
PAPER_NAME_START = '">'
PAPER_NAME_END = '</span>'
AUTHOR_BLOCK_START = '<div class="authors">'
AUTHOR_BLOCK_END = '<\div>'
ONE_AUTHOR_START = '<a href='
ONE_AUTHOR_END = '<\a>'
AUTHOR_START = 'data-author-name="'
AUTHOR_END = '"></span>'
AFFILIATION_START = 'data-affiliation="'
AFFILIATION_END = '|c|"></span>'
URL_START = '<div class="volumes">'
URL_END = '</div>'
# URL = 'http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=24214&punumber=9'
# Use below URL
URL = 'http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=24214&punumber=9'
#URL = 'http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=4456539&punumber=9&rowsPerPage=100'
# use below 2 urls to produce ruths network
URL1 = 'http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=6244883&rowsPerPage=100'
URL2 = 'http://ieeexplore.ieee.org/xpl/tocresult.jsp?isnumber=6564408&rowsPerPage=100'


G = zen.Graph()
visited = set()

# recursive function to find the nodes in all the proofwiki pages
def add_nodes(theurl):
	#print '\nAnalyzing: '+theurl
	try:
		f = urllib.urlopen(theurl)
		s = f.read()
		f.close()
	except:
		print '\t!!! Connection Issue'
		time.sleep(2)
		add_nodes(theurl)

	visited.add(theurl)

	contents = s.split(PAGE_START)[1] # keep only the part after the first occurance
	contents = contents.split(PAGE_END)[0]
	papers = contents.split(PAPER_START)[1:]
	for paper in papers:
		paper_content = paper.split(CONTENT_START)[1]
		paper_content = paper_content.split(CONTENT_END)[0]
		author_content = paper_content.split(AUTHOR_BLOCK_START)[1]
		author_content = author_content.split(AUTHOR_BLOCK_END)[0]
		if author_content:
			authors = author_content.split(ONE_AUTHOR_START)[1:]
			author_list = []
			affiliation_list = []
			for author in authors:	
				aff_data = author			
				author_name = author.split(AUTHOR_START)[1]
				author_name = author_name.split(AUTHOR_END)[0]
				#author_name = re.sub('[^A-Za-z0-9., ]+', '', author_name)
				soup = BeautifulSoup(aff_data, "html.parser")
				inputTags = []
				inputTags = soup.findAll("span")
				for tag in inputTags:
					if tag.has_attr('data-affiliation'):
 						affiliation = tag["data-affiliation"]
 						affiliation = affiliation.strip()
 						affiliation = affiliation.strip("|c|")
 						affiliation = re.sub('[&]', '', affiliation)
 						affiliation = re.sub('[""]', '', affiliation)
					if not affiliation:
						affiliation = "Independent Consultant"
				# add node only if author name is not empty
				if author_name:
					# author_name = author_name + " Affiliated to " + affiliation
					author_name = re.sub('[^A-Za-z0-9 ]+', '', author_name)
					author_list.append(author_name)
					affiliation_list.append(affiliation)
			# add nodes and edges only if author_list is not empty
			if len(author_list) >= 1:
				title_content = paper_content.split(TITLE_START)[1]
				title_content = title_content.split(TITLE_END)[0]
				soup = BeautifulSoup(title_content, "html.parser")				
				tag=soup.span
				if tag:
					paper_name = tag.string					
					# strip away all unwanted characters from author names
					author_list = [each_author.strip() for each_author in author_list]
					for each_author, affiliated_univ in zip(author_list, affiliation_list):
						if each_author not in G:
							G.add_node(each_author, data=affiliated_univ)
					# ADD EDGE ONLY IF THERE ARE MORE THAN 1 AUTHOR IN THE PAPER
					if len(author_list) > 1:
						# list all possible combinations of authors of a paper
						author_comb = list(itertools.combinations(author_list,2))
						# establish edge for each pair of combination
						for author_pair in author_comb:
							first_author = author_pair[0]
							second_author = author_pair[1]
							if G.has_edge(first_author,second_author):
								eidx = G.edge_idx(first_author,second_author)
								G.set_weight_(eidx,G.weight_(eidx)+1)
							else:
								G.add_edge(first_author,second_author,weight=1)	
	return G

def get_urls(theurl):
	print '\nAnalyzing: '+theurl
	print '\nGetting the URLs for past issues...'
	try:
		f = urllib.urlopen(theurl)
		s = f.read()
		f.close()
	except:
		print '\t!!! Connection Issue'
		time.sleep(2)
		add_nodes(theurl)
	url_list = []
	url_contents = s.split(URL_START)[1] # keep only the part after the first occurance
	url_contents = url_contents.split(URL_END)[0]
	lis = url_contents.split('<li>')[1:]
	for li in lis:
		url = li.split('<a href="')[1]
		url = url.split('"')[0]
		post_addition = "&rowsPerPage=100"
		url += post_addition
		pre_addition = "http://ieeexplore.ieee.org"
		url = pre_addition + url
		url_list.append(url)
	return url_list

# test cases
url_list = []
url_list = get_urls(URL)
#url_list = url_list[:3]
# url_list.append(URL1)
# url_list.append(URL2)
print 'Totally %i volumes were published from beginning.' %len(url_list)
i = 1
for each_url in url_list:
	print '\n %i' %i
	G = add_nodes(each_url)
	i = i + 1
	time.sleep(1)
# put to sleep for 2 sec between iteration
G_Nodes = G.num_nodes
G_Edges = G.num_edges
#Printing the number of nodes and edges in the network (G) given in Fig 6.1a
print '#Nodes: %i, #Edges: %i' % (G_Nodes,G_Edges)

# Get rid of weird characters
for nidx in range(0,G_Nodes): 
	nobj = unicode( G.node_object(nidx), 'utf-8' )
	nobj = unicodedata.normalize('NFKD', nobj).encode('ascii','ignore')
	G.set_node_object_(nidx,nobj)
				
zen.io.gml.write(G,'author_nw.gml')
