from playwright.sync_api import sync_playwright

def test_ui():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1920, "height": 1080})
        page.goto("http://localhost:8501")

        # Wait for app to load completely
        page.wait_for_selector(".stApp", timeout=60000)

        # Click search to reveal content
        page.evaluate('window.scrollTo(0, document.body.scrollHeight)')

        # Take a screenshot to verify UI components visually
        page.screenshot(path="ui_test.png")
        print("UI test complete, screenshot saved.")

        browser.close()

if __name__ == "__main__":
    test_ui()