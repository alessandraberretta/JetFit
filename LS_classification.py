import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from time import sleep
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.ui import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager


path_dir = '/Users/alessandraberretta/2014_rebin_nf_results/'
summary_list = [path_dir +
                file for file in os.listdir(path_dir) if file.startswith('summary_')]

for file in summary_list:
    GRB_name = file.split('_')[4] + '_' + file.split('_')[5]
    if GRB_name == 'GRB_140903A':
        df = pd.read_csv(file, sep='\t', index_col=0)
        df2 = pd.DataFrame([['Class_T90', 'S']], columns=[
                           'Parameters', 'Values'])
        df3 = df.append(df2, sort=False, ignore_index=True)
        df3.to_csv('new_sum_' + GRB_name + '.csv',  sep='\t')
    else:
        df = pd.read_csv(file, sep='\t', index_col=0)
        df2 = pd.DataFrame([['Class_T90', 'L']], columns=[
                           'Parameters', 'Values'])
        df3 = df.append(df2, sort=False, ignore_index=True)
        df3.to_csv('new_sum_' + GRB_name + '.csv',  sep='\t')

'''
def waitFor(driver: "WebDriver", element: tuple, timeout: int = 10):
    try:
        element_present = expected_conditions.presence_of_element_located(
            element)
        WebDriverWait(driver, timeout).until(element_present)
    except TimeoutException:
        print("Timed out waiting for page to load")
    finally:
        print("Wait done...")


def getGRBdf(driver: "WebDriver") -> pd.core.frame.DataFrame:
    driver.get("https://www.swift.ac.uk/xrt_live_cat/")
    waitFor(driver, (By.ID, "indexTable"))
    dfs = pd.read_html(driver.page_source)
    df = dfs[0]
    return df[df['Redshift'] >= 0]


def getGRBList(driver: "WebDriver") -> list:
    driver.get("https://www.swift.ac.uk/xrt_live_cat/")
    waitFor(driver, (By.ID, "indexTable"))
    html = driver.page_source
    soup = BeautifulSoup(html, "html.parser")
    return soup.select("#indexTable tbody tr")


def get_data(driver: "WebDriver") -> str:
    html = driver.page_source
    soup = BeautifulSoup(html, "html.parser")
    return soup.select("pre")[0].get_text()


def hasKnownRedshift(grb_sr: pd.Series, grb: str = "") -> bool:
    for elm in grb_sr.str.contains(pat=grb):
        if elm:
            return True
    return False


def get_single_T90(driver: "WebDriver", grbdf: pd.core.frame.DataFrame, grb_sr: pd.Series, grb: str = ""):
    out_filename = grb
    T90 = ""
    if hasKnownRedshift(grb_sr, grb):

        grb_list = getGRBList(driver)

        for tr in grb_list:
            if tr.get_text().find(f"GRB {grb}") != -1:
                link = driver.find_element_by_link_text(tr.a.get_text())
                link.click()
                sleep(0.5)
                waitFor(driver, (By.ID, "related"))
                link = driver.find_element_by_link_text("BAT analysis")
                link.click()
                sleep(0.5)
                try:
                    BAT_report = driver.find_element_by_xpath(
                        '/html/body/pre[1]').text
                    T90 = BAT_report[BAT_report.find(
                        'T9') + 5:BAT_report.find('Measured from')]
                except NoSuchElementException:
                    print(f"Cannot get GRB {grb} T90")
    return out_filename, T90


def get_annual_T90(driver: "WebDriver", grbdf: pd.core.frame.DataFrame, grb_sr: pd.Series, year: str = ""):

    filenames = []
    T90 = []

    for tr in getGRBList(driver):
        if f"GRB {year}" in tr.get_text():
            tmp_filename, tmp_T90 = get_single_T90(driver, grbdf, grb_sr, str(
                tr.get_text()[tr.get_text().find(' ')+1: tr.get_text().find('h') - 2]))
            if tmp_filename and tmp_T90:
                filenames.append(tmp_filename)
                T90.append(tmp_T90)
    return filenames, T90


def main():

    chrome_options = Options()

    # chrome_options.add_argument("--headless")

    driver = webdriver.Chrome(
        ChromeDriverManager().install(), options=chrome_options)

    grbdf = getGRBdf(driver)
    grb_sr = pd.Series(grbdf['GRB'])

    filenames, T90 = get_annual_T90(driver, grbdf, grb_sr, year='13')
    print(filenames)
    print(T90)

    for idx, elm in enumerate(T90):
        if float(elm.split(' ')[0]) <= 4:
            print(filenames[idx])


if __name__ == "__main__":
    main()
'''
