from playwright.sync_api import sync_playwright

def test_ui():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto("http://localhost:8501")
        page.wait_for_selector("text=Run Vina Docking")

        # Click the first matching button
        page.click("button:has-text('Run Vina Docking')")

        # Wait for completion message
        page.wait_for_selector("text=Docking complete!", timeout=60000)

        # Check that the dataframe is present
        page.wait_for_selector("[data-testid='stDataFrame']")

        # Check that st.text_area is not there
        assert page.locator("text=Vina Output").count() == 0

        page.screenshot(path="screenshot.png")
        browser.close()

if __name__ == "__main__":
    test_ui()
