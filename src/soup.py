from bs4 import BeautifulSoup
import urllib


def get_id(html, el, el_id):
    return html.find(el, id=el_id)


def no_result(j):
    return j['iTotalRecords'] == 0


def parse(v):
    return BeautifulSoup(v, "html.parser")


def j_data(j):
    return j['aaData']


def get_param(href, param):
    url = urllib.parse.urlparse(href)
    params = urllib.parse.parse_qs(url.query)
    if param in params:
        return params[param]
    return None
