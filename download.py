from enum import Enum
from time import sleep
import pandas as pd
import typer
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.ui import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager


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


'''
class Command(str, Enum):
    get_grb = "getGRB"
    download_grb = "download"
'''


def download_single_grb(driver: "WebDriver", grbdf: pd.core.frame.DataFrame, grb_sr: pd.Series, grb: str = "", ergcm2s: bool = False):
    print(f"Downloading data for GRB {grb}")
    out_filename = ""
    redshift = -999

    if hasKnownRedshift(grb_sr, grb):

        redshift = (grbdf['Redshift'][grbdf['GRB']
                                      == "GRB " + grb]).values[0]
        grb_list = getGRBList(driver)
        for tr in grb_list:
            if tr.get_text().find(f"GRB {grb}") != -1:
                link = driver.find_element_by_link_text(tr.a.get_text())
                link.click()
                sleep(0.5)
                waitFor(driver, (By.ID, "related"))
                if ergcm2s:
                    link = driver.find_element_by_link_text("Light curve")
                    link.click()
                    sleep(0.5)
                    try:
                        link = driver.find_element_by_id(
                            "flux_makeDownload")
                        link.click()
                        sleep(0.5)
                        waitFor(driver, (By.TAG_NAME, "pre"))
                        data = get_data(driver)
                        out_filename = f"GRB_{grb}_{redshift}.data"
                        with open(out_filename, "w") as out_file:
                            out_file.write(data)
                        print(f"Data written in {out_filename}")
                    except NoSuchElementException:
                        print(f"Cannot download GRB {grb} data...")
                else:

                    link = driver.find_element_by_link_text("Burst Analyser")
                    link.click()
                    sleep(0.5)
                    waitFor(driver, (By.ID, "XRT"))
                    try:
                        link = driver.find_element_by_id(
                            "xrt_DENSITY_makeDownload")
                        link.click()
                        sleep(0.5)
                        waitFor(driver, (By.TAG_NAME, "pre"))
                        data = get_data(driver)
                        out_filename = f"GRB_{grb}_{redshift}.data"
                        with open(out_filename, "w") as out_file:
                            out_file.write(data)
                        print(f"Data written in {out_filename}")
                    except NoSuchElementException:
                        print(f"Cannot download GRB {grb} data...")

    return out_filename, redshift


def download_year_grb(driver: "WebDriver", grbdf: pd.core.frame.DataFrame, grb_sr: pd.Series, year: str = "", ergcm2s: bool = False):
    print(f"Downloading 20{year} GRB data")
    filenames = []
    redshifts = []

    for tr in getGRBList(driver):
        if f"GRB {year}" in tr.get_text():
            tmp_filename, tmp_redshift = download_single_grb(driver, grbdf, grb_sr, str(
                tr.get_text()[tr.get_text().find(' ')+1: tr.get_text().find('h') - 2]), ergcm2s)
            if tmp_filename and tmp_redshift != -999:
                filenames.append(tmp_filename)
                redshifts.append(tmp_redshift)

    return filenames, redshifts


def download_data(grb: str = "", year: str = "", headless: bool = False, ergcm2s: bool = False):
    chrome_options = Options()

    if headless:
        chrome_options.add_argument("--headless")

    driver = webdriver.Chrome(
        ChromeDriverManager().install(), options=chrome_options)

    '''
    if command == "getGRB":
        print("\nShowing GRB list with known redshift...\n")
        print(getGRBdf(driver))
    elif command == "download":
    '''

    grbdf = getGRBdf(driver)
    grb_sr = pd.Series(grbdf['GRB'])

    if grb != "":
        # Search a specific GRB
        filename, redshift = download_single_grb(
            driver, grbdf, grb_sr, grb, ergcm2s)
    elif year != -1:
        # Search all GRBs in a specific year
        filename, redshift = download_year_grb(
            driver, grbdf, grb_sr, year, ergcm2s)

    driver.close()

    return filename, redshift


'''
if __name__ == "__main__":
    typer.run(main)
'''
